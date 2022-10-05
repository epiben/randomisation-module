library(R6)
library(glue)
library(purrr)
library(lubridate)
library(dplyr)

PermutedBlockList <- R6Class("PermutedBlockList",
  public = list(
    queue = list(),
    max_allocs = NULL,
    n_randomised = 0,
    arms = NULL,
    block_sizes = NULL,
    enque_n_blocks = NULL,
    seed = NULL,

    initialize = function(arms = c("A", "B"), block_sizes = c(4, 6, 8),
                          enque_n_blocks = 3, max_allocs = Inf, seed = 4131) {
      self$max_allocs <- max_allocs
      self$arms <- arms
      self$block_sizes <- block_sizes
      self$enque_n_blocks <- enque_n_blocks
      self$seed <- seed

      private$n_arms <- length(arms)

      private$check_block_params()
      private$enque()
      self$show_status()
    },

    show_status = function() {
      cat(glue(
        "Sizes: {paste(self$block_sizes, collapse = ', ')}",
        "No. block to add: {self$enque_n_blocks}",
        "Allocations remaining: {length(self$queue)}",
        .sep = " -- "
      ))
    },
    allocate_next = function() {
      # NB: Inefficient implementation, but fine for prototyping with short lists
      if (self$n_randomised >= self$max_allocs) {
        cat("Will randomise no more patients, max_allocs reached")
        return(NULL)
      }
      if (length(self$queue) == 0) private$enque()

      val <- self$queue[[1]]
      self$queue <- self$queue[-1]
      self$n_randomised <- self$n_randomised + 1
      return(val$arm)
    },
    show_tidy_list = function() {
      bind_rows(self$queue)
    }
  ),

  private = list(
    n_arms = NULL,
    allocation_id = 1,
    block_id = 1,

    check_block_params = function() {
      if (any(self$block_sizes %% length(self$arms) != 0)) {
        stop("All blocks sizes must be a multiple of the number of arms", call. = FALSE)
      }
    },
    set_seed = function() {
      seed_used <- self$seed
      set.seed(seed_used)
      self$seed <- seed_used + 1
      invisible(seed_used)
    },
    enque = function() {
      allocate_within_block <- function(block_size, arms) {
        allocs <- list()
        for (i in seq_len(block_size)) {
          seed_used <- private$set_seed() # sets a new seed and returns it
          weights <- block_size/private$n_arms - map_dbl(arms, ~ sum(. == allocs))
          allocs[[i]] <- list(
            allocation_id = private$allocation_id,
            block_id = private$block_id,
            arm = sample(arms, 1, prob = weights),
            seed_used = seed_used
          )
          private$allocation_id <- private$allocation_id + 1
        }
        private$block_id <- private$block_id + 1
        stopifnot(all(table(allocs$arm) == block_size/private$n_arms))
        return(allocs)
      }

      private$set_seed()
      blocks <- sample(self$block_sizes, self$enque_n_blocks, TRUE)
      self$queue <- flatten(map(blocks, allocate_within_block, arms = self$arms))
    }
  )
)

q1 <- PermutedBlockList$new()
q1$show_tidy_list()

for (i in seq_len(30)) {
  print(glue("{i}: {q1$allocate_next()}"))
}

q2 <- PermutedBlockList$new(arms = c("A", "B", "C"), block_sizes = c(6, 9, 12), seed = 1)
as.data.frame(q2$show_tidy_list())


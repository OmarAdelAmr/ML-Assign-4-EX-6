directory <- '~/Desktop/ML-Assign-4-EX-6'
setwd(directory)

training_set <- read.csv('Sequences_train.csv', sep = ',', header = FALSE, stringsAsFactors = FALSE)
testing_set <- read.csv('Sequences_test_unlabeled.csv', sep = ',', header = FALSE, stringsAsFactors = FALSE)

subsquence_length <- 3

get_all_subsequences <- function(k, input_string)
{
  result <- substring(input_string, 1:(nchar(input_string) - k + 1), k:nchar(input_string))
  
  return(result)
}


count_occurances <- function(subsquence, subsequences_list)
{
  return(sum(subsequences_list == subsquence))
}

matrix_size = nrow(training_set)
A = matrix(nrow = matrix_size, ncol = matrix_size, byrow = TRUE)

for (x in 1:matrix_size)
{
  # temp_dataset <- training_set[-c(x),]
  for (y in 1:matrix_size)
  {
      x_subsequences = get_all_subsequences(subsquence_length, training_set[x, 1])
      y_subsequences = get_all_subsequences(subsquence_length, training_set[y, 1])
      match_score <- 0
      for (sub in x_subsequences)
      {
        match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
      }
      A[x, y] <- match_score
  }
}

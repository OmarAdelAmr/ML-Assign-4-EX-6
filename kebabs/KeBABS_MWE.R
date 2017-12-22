## title: Example: Pre-Computed Kernel Matrices with SVM Optimization via KeBABS
## author: Ulrich Bodenhofer

## to run this file as a whole and produce a PDF report, enter the following:
## > install.packages("knitr") ## if not already installed
## > library(knitr)
## > stitch("KeBABS_MWE.R")

library(kebabs)
## if the package is not installed, enter the following:
## > source("https://bioconductor.org/biocLite.R")
## > biocLite("kebabs")

## load data
data(iris)

## transform into binary -1/+1 classification problem
iris$Species <- factor(ifelse(iris$Species == "setosa", +1, -1))

train <- sort(sample(1:nrow(iris), 0.6 * nrow(iris)))
test <- (1:nrow(iris))[-train]

## compute simple kernel matrix using linear kernel
Kmat <- tcrossprod(as.matrix(iris[train, 1:4]))

## normalize kernel matrix
selfSim <- sqrt(diag(Kmat))
Kmat <- sweep(Kmat, 1, selfSim, FUN="/") # scale rows
Kmat <- sweep(Kmat, 2, selfSim, FUN="/") # scale columns
Kmat <- as(Kmat, "KernelMatrix")

## train SVM
model <- kebabs:::svmd.default(Kmat, iris$Species[train], cost=1)

## compute "kernel matrix" of test vs. support vectors
## (model$index contains indices of support vectors)
KmatTest <- tcrossprod(as.matrix(iris[test, 1:4]),
                       as.matrix(iris[train[model$index], 1:4]))
selfSimTest <- sqrt(apply(as.matrix(iris[test, 1:4]), 1,
                          function(x) sum(x^2)))
KmatTest <- sweep(KmatTest, 1, selfSimTest, FUN="/") # scale rows
KmatTest <- sweep(KmatTest, 2, selfSim[model$index], FUN="/") # scale columns
KmatTest <- as(KmatTest, "KernelMatrix")

## make predictions for test samples
## (model$coefs  ... alpha_i * y_i for support vectors,
##  model$rho    ... -b,
##  model$levels ... labels used by SVM,
##  model$label  ... assignment of labels)
yTest <- factor(model$levels[ifelse(KmatTest %*% model$coefs - model$rho >= 0,
                                    model$label[1], model$label[2])])

## compute confusion table
table(y=iris$Species[test], g=yTest)

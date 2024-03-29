#install packages
install.packages("tidyverse")
install.packages("dplyr")
install.packages("MASS")
install.packages("Matrix")
install.packages("stats4")
install.packages("gmm")
install.packages("sandwich")
install.packages("tmvtnorm")
install.packages("crch")
install.packages("coda")
install.packages("MCMCpack")
library(tmvtnorm)
library(mvtnorm)
library(Matrix)
library(stats4)
library(gmm)
library(sandwich)
library("stats")

# data preprocessing steps--------------------------------------------------------------------
set.seed(100)
x1<-matrix(c(rep(1, 100)), ncol=1)
x2_x3 <- matrix(rnorm(200), nrow = 100, ncol = 2)
x <- cbind(x1, x2_x3)
eps<-matrix(rnorm(100),nrow =100, ncol=1)
beta<- matrix(c(1, 2, 0.5), ncol = 1)
y<-x%*%beta+eps




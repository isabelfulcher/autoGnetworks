library(readr)
library(devtools)
library(autognet)
library(tidyverse)
library(MASS)

args <- commandArgs(trailingOnly = TRUE)
samplesize <- args[1] 
type <- args[2] #coding or pl 
iter <- as.numeric(args[3])


#############################
#######   LOAD DATA   #######
#############################

#functions
source("01_functions.R")

#truth 
load(paste0("input/truth_",samplesize,".rda"))
load(paste0("input/truth_causal_estimands_",samplesize,".rda"))


#Gibbs results 
results <- readRDS(paste0("output/gibbs_",type,"_",samplesize,"_all.rds"))

#graph characteristics
file_name_adjmat=paste0("input/adj",as.character(samplesize))

file_name_maxindep=paste0("input/maxindep",as.character(samplesize))

adjmat <- as.matrix(read_delim(file_name_adjmat, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
maxindep <- t(as.matrix(read_delim(file_name_maxindep, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)))

#Turn adjacency matrix into list of neighbors for each individual (needed for sims)
N <- nrow(adjmat)
adjacency <- list(NA)
for (i in 1:N){
  adjacency[[i]] <- which(adjmat[i,]==1)
}

weights <- apply(adjmat,1,sum)

#############################
#######   ESTIMATES   #######
#############################

#Input values 
alpha <- results[[1]][iter,]
alpha.vcov <- results[[5]][[iter]]

beta <- results[[2]][iter,]
beta.vcov <- results[[6]][[iter]]

group_lengths <- c(1,1,1)
group_functions <- c(1,1,1)

pr_trt <- 0.7

## Make object for AGC pkg 
beta.new <- beta[c(1,2,4,6,8,10,3,5,7,9)] #reorder because of coding function output
alpha.new <- alpha[c(1:9,11,10,12,14,15,13)] #reorder because of coding function output

outlist.point <- list(t(as.matrix(alpha.new)),t(as.matrix(beta.new)),
                NA,NA,NA,
                group_lengths,group_functions,adjmat)
names(outlist.point) <- c("alpha", "beta", "NA", "NA", "NA", "group_lengths", "group_functions", "adjmat")
class(outlist.point) <- append(class(outlist.point),"agcParamClass")

## Run AGC package 
R <- 1000
burnin_R <- 200

point.estimate <- agcEffect(outlist.point, burnin = 0, thin = 1, treatment_allocation = pr_trt, subset = 0,
                            R = R, burnin_R = burnin_R, burnin_cov = 0, average = TRUE, index_override = 0,
                            return_effects = 0)

## Save values 
results.point.estimate <- data.frame(average = point.estimate[1],
                                     overall = point.estimate[1] - point.estimate[2],
                                     direct = point.estimate[3] -  point.estimate[4],
                                     spillover = point.estimate[4] -  point.estimate[2])

#############################
#######   VARIANCE    #######
#############################

if(type == "coding"){

  B <- 200
  set.seed(10*iter) 
  draw.l <- MASS::mvrnorm(B,alpha,alpha.vcov)
  draw.y <- MASS::mvrnorm(B,beta,beta.vcov)
  draw.y.new <- draw.y[,c(1,2,4,6,8,10,3,5,7,9)] #reorder because of coding function output
  draw.l.new <- draw.l[,c(1:9,11,10,12,14,15,13)] #reorder because of coding function output
  
  outlist <- list(draw.l.new,draw.y.new,
                  NA,NA,NA,
                  group_lengths,group_functions,adjmat)
  names(outlist) <- c("alpha", "beta", "NA", "NA", "NA", "group_lengths", "group_functions", "adjmat")
  class(outlist) <- append(class(outlist),"agcParamClass")
  
  
  estimates.boot <- agcEffect(outlist, burnin = 0, thin = 1, treatment_allocation = pr_trt, subset = 0,
                              R = 50, burnin_R = 0, burnin_cov = 10, average = TRUE, index_override = 0,return_effects = 0)
  ## Save values
  results.boot <- data.frame(average = estimates.boot[,1],
                             overall = estimates.boot[,1] - estimates.boot[,2], 
                             direct = estimates.boot[,3] -  estimates.boot[,4],
                             spillover = estimates.boot[,4] -  estimates.boot[,2])
  
  results.median.boot <- apply(results.boot,2,median)
  results.var.boot <- apply(results.boot,2,var)
  
  ## Coverage 
  average.truth <- apply(check.truth,2,mean)[1]
  coverage.boot <- data.frame(average = as.numeric( (quantile(results.boot$average,.025) <= average.truth) & (quantile(results.boot$average,.975) >= average.truth) ),
                              overall = as.numeric( (quantile(results.boot$overall,.025) <= truth.final[1]) & (quantile(results.boot$overall,.975) >= truth.final[1]) ),
                              direct = as.numeric( (quantile(results.boot$direct,.025) <= truth.final[2]) & (quantile(results.boot$direct,.975) >= truth.final[2]) ),
                              spillover = as.numeric( (quantile(results.boot$spillover,.025) <= truth.final[3]) & (quantile(results.boot$spillover,.975) >= truth.final[3]) ))
  
  coverage.wald <- function(x,var,truth) { as.numeric((x + qnorm(.025)*sqrt(var) <= truth) & (x + qnorm(.975)*sqrt(var) >= truth)) }
  coverage.boot.wald <- data.frame(average = coverage.wald(results.point.estimate[1],results.var.boot[1],average.truth),
                                   overall = coverage.wald(results.point.estimate[2],results.var.boot[2],truth.final[1]),
                                   direct = coverage.wald(results.point.estimate[3],results.var.boot[3],truth.final[2]),
                                   spillover = coverage.wald(results.point.estimate[4],results.var.boot[4],truth.final[3]))
}

#############################
#########   SAVE   ##########
#############################

if (type == "coding"){
  results <- list(results.point.estimate,results.median.boot,results.var.boot,coverage.boot,coverage.boot.wald)
} else {
  results <- list(results.point.estimate)
}

saveRDS(results, file = paste0("output/estimates_",type,"_",samplesize,"_iter",iter,".rds"))

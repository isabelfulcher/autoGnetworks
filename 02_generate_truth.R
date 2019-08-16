library(readr)
library(autognet)

args <- commandArgs(trailingOnly = TRUE)
samplesize <- args[1]

#############################
#######   LOAD DATA   #######
#############################

load(paste0("input/truth_",as.character(samplesize),".rda"))

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

#####################################
#######   GENERATE DATASETS   #######
#####################################
#STEP 0. REORDER TRUTH VECTOR FOR IMPLEMENTATION INTO COVER GENERATION
beta.truth.2 <- beta.truth[c(1,2,4,6,8,10,3,5,7,9)]

#STEP 1. NECESSARY INPUTS 
pr_trt <- .7
R <- 3000
burnin_R <- 1000
group_lengths <- rep(1,3)
group_functions <- rep(1,3)

#STEP 2. CREATE RESULT OBJECT FOR AGC FUNCTION 
alpha.truth <- rbind(c(tau.truth,rho.truth,nu.truth),c(tau.truth,rho.truth,nu.truth),c(tau.truth,rho.truth,nu.truth))
beta.truth <- rbind(beta.truth.2,beta.truth.2,beta.truth.2)
outlist <- list(alpha.truth,beta.truth,
                NA,NA,NA,
                group_lengths,group_functions,adjmat)
names(outlist) <- c("alpha", "beta", "NA", "NA", "NA", "group_lengths", "group_functions", "adjmat")
class(outlist) <- append(class(outlist),"agcParamClass")

#STEP 3. TRUTH ESTIMATES 

check.truth <- agcEffect(outlist, burnin = 0, thin = 1, treatment_allocation = pr_trt, subset = 0,
                         R = R, burnin_R = burnin_R, burnin_cov = 0, average = TRUE, index_override = 0, 
                         return_effects = 0)

#STEP 4. CALCULATE CAUSAL EFFECTS 
overall <- check.truth[,1] - check.truth[,2]
direct <- check.truth[,3] - check.truth[,4]
spillover <- check.truth[,4] - check.truth[,2]
truth <- cbind(overall,direct,spillover)
truth.final <- apply(truth,2,mean)

##################################
############   SAVE   ############
##################################

save(truth.final,check.truth,
     file = paste0("input/truth_causal_estimands_",samplesize,".rda"))

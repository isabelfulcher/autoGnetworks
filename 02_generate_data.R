library(readr)

args <- commandArgs(trailingOnly = TRUE)
samplesize <- args[1]

source("01_functions.R")

#############################
##  TRUE PARAMETER VALUES  ##
#############################
# NOTE: Follows the structure laid out in the AGC paper 

#covariate model
tau.truth <- c(-1,0.5,-0.5) #intercepts for L1, L2, L3 
rho.truth <- c(.1,.2,.1) #how other covariates affect current covariate (1-2, 1-3, 2-3)
nu.truth <- c(.1,.1,.1) #NOTE: technically c(.1,0,0,.1,0,0,.1,0,0) (1-1, 1-2, 1-3, 2-2, 2-1, 2-3, 3-3, 3-1, 3-2) 

#treatment model
gamma.truth <- c(-1,0.2,0.3,0.2,0.3,-0.2,-0.3,0.7)

#outcome model
beta.truth <- c(-0.3,-1.0,-0.9,0.1,0.2,0.1,0.2,-0.1,-0.2,0.3)

#############################
#######   LOAD DATA   #######
#############################

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
gamma.truth.2 <- gamma.truth[c(1,2,4,6,8,3,5,7)]
beta.truth.2 <- beta.truth[c(1,2,4,6,8,10,3,5,7,9)]

#STEP 1. 
S <- 1000 #number of samples
burnin <- 1000 
thin <- 5
R <- S*thin + burnin
group_lengths <- rep(1,3)
group_functions <- rep(1,3)

#STEP 2. GENERATE COVARIATE ARRAY 
cov.array <- gen.var.covariates2(tau.truth, rho.truth, nu.truth, 
                                 N, R, adjacency, weights, group_lengths, group_functions)

#STEP 3. GENERATE TREATMENT
trt.gen <- gen.var.treatment2(cov.array,gamma.truth.2,N,R,adjacency,weights)

#STEP 4. GENERATE OUTCOME
outcome.gen <- gen.var.outcome2(trt.gen,cov.array,beta.truth.2,N,R,adjacency,weights)

#STEP 5. CONSOLIDATE INTO DATASETS
cov.final <- cov.array[,,(burnin+1):R][,,seq(1,S*thin,by=thin)]
trt.final <- trt.gen[,(burnin+1):R][,seq(1,S*thin,by=thin)]
outcome.final <- outcome.gen[,(burnin+1):R][,seq(1,S*thin,by=thin)]

data.gen <- lapply(1:S,function(s) cbind(outcome.final[,s],trt.final[,s],cov.final[,,s]))


##################################
############   SAVE   ############
##################################

save(data.gen, tau.truth, rho.truth, nu.truth, gamma.truth, beta.truth, 
     file = paste0("input/truth_",samplesize,".rda"))


##misspecfied graphs have same truth -- obviates downstream issues
if (samplesize == "800") {
  
  save(data.gen, tau.truth, rho.truth, nu.truth, gamma.truth, beta.truth, 
       file = paste0("input/truth_800miss1.rda"))
  
  save(data.gen, tau.truth, rho.truth, nu.truth, gamma.truth, beta.truth, 
       file = paste0("input/truth_800miss2.rda"))
  
}

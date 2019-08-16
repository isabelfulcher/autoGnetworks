###########################################
######    DEFINE FUNTIONS FOR SIMS    #####
###########################################

#################################
#### GENERATE MARKOV COVER ######
#################################

gen.var.covariates2 <- function(tau, rho, nu, N, R,
                                adjacency, weights, group_lengths, group_functions){
  
  #x is a dummy variable so i can apply over it
  J <- length(tau)
  
  cov.mat <- sapply(1:J,function(x) rbinom(N,1,runif(1,.1,.9)))
  cov_save <- array(NA,c(N,J,R))
  
  # Make symmetrical matrix of the rho values
  rho_mat <- matrix(0, nrow = J, ncol = J)    
  rho_mat[lower.tri(rho_mat, diag=FALSE)] <- rho; rho_mat <- rho_mat + t(rho_mat)
  
  # Number of iterations
  for (r in 1:R){
    
    # Number of people
    for (i in 1:N){
      
      # Index covariate 
      j <- 1
      
      # Number of groups (of covariants)
      for (group_index in 1:length(group_lengths)){
        
        # Values associated with each group
        group_length <- group_lengths[group_index]
        group_function <- group_functions[group_index]
        
        # group_length is the number of binarized covariates 
        # group_function is a numeric dictionary key that utilizes a value from the covariate_process function
        # if group_length is > 1, meaning that we have multiple binarized covariates associated with a specific
        # entity, then it is definitely a multinomial
        if(group_length > 1){
          
          # Multinomial case
          prob_vec <- sapply(0:(group_length-1), function(m){
            j_prime <- j + m
            exp(tau[j_prime] + sum(rho_mat[,j_prime]*cov.mat[i,]) + nu[j_prime]*sum(cov.mat[adjacency[[i]],j_prime]/weights[i]))
          })
          
          # for the rmultinom call, have to append a 1 and remove the last value; update several values
          cov.mat[i, j + (0:(group_length-1))] <- (rmultinom(1,1,c(prob_vec,1))[,1])[-1*(group_length+1)]
          
        } else if(group_function == 1){
          
          # Logistic / binary case
          prob_Lj <- plogis(tau[j] + sum(rho_mat[,j]*cov.mat[i,]) + nu[j]*sum(cov.mat[adjacency[[i]],j]/weights[i]))
          cov.mat[i,j] <- rbinom(1,1,prob_Lj)
        } # add in normal here once form is decide
        
        j <- j + group_length
        
      }
    }
    cov_save[,,r] <- cov.mat
  }
  return(cov_save)
}


gen.var.treatment2 <- function(cov,gamma,N,R,adjacency,weights){
  # gamma is vector of parameters for tretamtnet model 
  # cov is a Nxp matrix of covariate values
  # N is number of individuals in network THAT HAVE AT LEAST ONE NEIGHBOR
  # R is number of iterations to run through
  # adjacency is a list of length N with neighbor IDs for each individual
  # start is starting point for chain (start at observed value of the treatment) 
  # weights is a vector of length N with the number of neighbors for each individual
  
  cov_all <- cov
  trt_save <- matrix(NA,N,R)
  
  vec <- rbinom(N,1,runif(1,.1,.9))
  prob_y <- NA
  
  for (r in 1:R){
    cov <- cov_all[,,r]
    for (i in 1:N){
      prob_outcome <- plogis(gamma[1] + # intercept term
                               gamma[2:(1+ncol(cov))]%*%cov[i,] + # individual covariate term(s)
                               gamma[(2+ncol(cov))]*sum(vec[adjacency[[i]]]/weights[i]) + # neighbor treatment
                               sum(gamma[(2+ncol(cov)+(1:ncol(cov)))]*colSums(cov[adjacency[[i]],(1:ncol(cov)),drop = FALSE]/weights[i])) # neighbor covariate
      )
      vec[i] <- rbinom(1,1,prob_outcome)
    }
    trt_save[,r] <- vec
  }
  return(trt_save)
}


gen.var.outcome2 <- function(trt,cov,beta,N,R,adjacency,weights){
  # beta is vector of parameters for outcome model (beta.p in 04a) 
  # trt is Nx1 matrix of treatment values (trt.i in 04a)
  # cov is a Nxp matrix of covariate values (cov.i in 04a) 
  # N is number of individuals in network THAT HAVE AT LEAST ONE NEIGHBOR
  # R is number of iterations to run through
  # adjacency is a list of length N with neighbor IDs for each individual
  # start is starting point for chain (outcome.i in 04a)
  # weights is a vector of length N with the number of neighbors for each individual
  
  #pastis the treatment and covariate chains generated previously (list of length 2) 
  trt_all <- trt
  cov_all <- cov
  outcome_save <- matrix(NA,N,R)
  
  vec <- rbinom(N,1,runif(1,.1,.9)) 
  prob_y <- NA
  
  for (r in 1:R){
    trt <- trt_all[,r]
    cov <- cov_all[,,r]
    
    for (i in 1:N){
      prob_outcome <- plogis(beta[1] + # intercept term
                               beta[2]*trt[i] + # individual treatment term
                               beta[3:(2+ncol(cov))]%*%cov[i,] + # individual covariate term(s)
                               beta[(3+ncol(cov))]*sum(vec[adjacency[[i]]]/weights[i]) + # neighbor outcome 
                               beta[(4+ncol(cov))]*sum(trt[adjacency[[i]]]/weights[i]) + # neighbor treatment
                               sum(beta[(4+ncol(cov)+(1:ncol(cov)))]*colSums(cov[adjacency[[i]],(1:ncol(cov)),drop = FALSE]/weights[i])) # neighbor covariate
      )
      
      # Isabel's sapply for the neighbor covariate line above
      # sum(sapply(1:ncol(cov),function(x){beta[(4+ncol(cov)+x)]*colSums(cov[adjacency[[i]],x]/weights[i])}))
      
      vec[i] <- rbinom(1,1,prob_outcome)
    }
    outcome_save[,r] <- vec
  }
  return(outcome_save)
}


#################################
############ CODING #############
#################################

## SCORE FUNCTION (3 COVARIATES) ##
U <- function(estimates,covariate,covariate.n){
  
  tau <- estimates[1:3]
  rho <- estimates[4:6]
  nu <- estimates[7:15]
  
  mean.L1 <- plogis(cbind(1,covariate[,2],covariate[,3],covariate.n[,1],covariate.n[,2],covariate.n[,3])%*%c(tau[1],rho[1],rho[2],nu[1],nu[2],nu[3]))
  mean.L2 <- plogis(cbind(1,covariate[,1],covariate[,3],covariate.n[,1],covariate.n[,2],covariate.n[,3])%*%c(tau[2],rho[1],rho[3],nu[5],nu[4],nu[6]))
  mean.L3 <- plogis(cbind(1,covariate[,1],covariate[,2],covariate.n[,1],covariate.n[,2],covariate.n[,3])%*%c(tau[3],rho[2],rho[3],nu[8],nu[9],nu[7]))
  
  
  score <- cbind( cbind(1,covariate[,2],covariate[,3],covariate.n[,1],covariate.n[,2],covariate.n[,3])*c(covariate[,1] - mean.L1),
                  cbind(1,covariate[,1],covariate[,3],covariate.n[,1],covariate.n[,2],covariate.n[,3])*c(covariate[,2] - mean.L2),
                  cbind(1,covariate[,1],covariate[,2],covariate.n[,1],covariate.n[,2],covariate.n[,3])*c(covariate[,3] - mean.L3))
  
  score.new <- score[,c(1,7,13,2,8,3,14,9,15,4,5,6,11,10,12,18,16,17)]
  
  score.new2 <- cbind(score.new[,1:3],score.new[,4]+score.new[,5],score.new[,6]+score.new[,7],score.new[,8]+score.new[,9],score.new[,10:18])
  
  return(score.new2)
}

## PSEUDOLIKELIHOOD (3 COVARIATES) ##
cov.pl <- function(parameters,covariate,covariate.n){
  
  tau <- parameters[1:3]
  rho <- parameters[4:6]
  nu <- parameters[7:15]
  
  quasilik <- sum( log((plogis(tau[1] + rho[1]*covariate[,2] + rho[2]*covariate[,3] + nu[1]*covariate.n[,1] + nu[2]*covariate.n[,2] + nu[3]*covariate.n[,3])^covariate[,1])
                       *(1-plogis(tau[1] + rho[1]*covariate[,2] + rho[2]*covariate[,3] + nu[1]*covariate.n[,1] + nu[2]*covariate.n[,2] + nu[3]*covariate.n[,3]))^(1-covariate[,1]))
                   + log((plogis(tau[2] + rho[1]*covariate[,1] + rho[3]*covariate[,3] + nu[5]*covariate.n[,1] + nu[4]*covariate.n[,2] + nu[6]*covariate.n[,3])^covariate[,2])
                         *(1-plogis(tau[2] + rho[1]*covariate[,1] + rho[3]*covariate[,3] + nu[5]*covariate.n[,1] + nu[4]*covariate.n[,2] + nu[6]*covariate.n[,3]))^(1-covariate[,2]))
                   +  log((plogis(tau[3] + rho[2]*covariate[,1] + rho[3]*covariate[,2] + nu[8]*covariate.n[,1] + nu[9]*covariate.n[,2] + nu[7]*covariate.n[,3])^covariate[,3])
                          *(1-plogis(tau[3] + rho[2]*covariate[,1] + rho[3]*covariate[,2] + nu[8]*covariate.n[,1] + nu[9]*covariate.n[,2] + nu[7]*covariate.n[,3]))^(1-covariate[,3]))
  )
  
  return(-1*quasilik)
}

## CONVERGENCE ISSUES (FOR SMALL MAXINDEP SETS IN SIMS) 

withWarnings  <- function (expr) { 
  warnings <- character() 
  retval <- withCallingHandlers(expr, warning = function(ex) { 
    warnings <<- c(warnings, conditionMessage(ex)) 
    invokeRestart("muffleWarning") 
  }) 
  list(Value = retval, Warnings = warnings) 
} 


####################################
############ BOOTSTRAP #############
####################################

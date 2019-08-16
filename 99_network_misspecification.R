library(readr)

samplesize <- "800" #removing from dense network

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
#######   CONVOLUTE NETWORK   #######
#####################################

## MISSING EDGES ONLY
set.seed(20)
ids <- sample(1:N,200)
miss.list <- sapply(ids,function(x) { sample(adjacency[[x]],sample(1:5,1)) })

adjmat.miss1 <- adjmat
for (i in 1:length(ids)){
  
  adjmat.miss1[miss.list[[i]],ids[i]] <- 0 
  adjmat.miss1[ids[i],miss.list[[i]]] <- 0 
  
}

weights.miss1 <- apply(adjmat.miss1,1,sum)


## MISSING NODES (WHICH IMPLIES MISSING EDGES) 
set.seed(10)
ind <- sample(1:N,75,replace=FALSE)
adjmat.miss2 <- adjmat[-ind,-ind]

weights.miss2 <- apply(adjmat.miss2,1,sum)


###########################################
#######   SAVE ADJACENCY MATRICES   #######
###########################################

write.table(adjmat.miss1, file = "revisions/input/adj800miss1.txt",col.names = F, row.names = F) 
write.table(adjmat.miss2, file = "revisions/input/adj800miss2.txt",col.names = F, row.names = F) 

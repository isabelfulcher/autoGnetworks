library(readr)
library(BuenColors)
library(tidyverse)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
samplesize <- args[1]

#############################
#######   LOAD DATA   #######
#############################

#truth
load(paste0("input/truth_",samplesize,".rda"))

#coding
results.coding <- readRDS(paste0("output/gibbs_coding_",samplesize,"_all.rds"))

#pl
results.pl <- readRDS(paste0("output/gibbs_pl_",samplesize,"_all.rds"))

#######################################
#######  CONSOLIDATE RESULTS    #######
#######################################

## ALPHA ##
nu.truth.full <- c(nu.truth[1],0,0,nu.truth[2],0,0,nu.truth[3],0,0)
alpha.truth <- c(tau.truth,rho.truth,nu.truth.full)

# coding
alpha.convergence <- results.coding[[9]]
alpha.values.coding <- results.coding[[1]][alpha.convergence==0,]
alpha.est.coding <- apply(results.coding[[1]][alpha.convergence==0,],2,mean,na.rm=TRUE)
alpha.mcvar.coding <- apply(results.coding[[1]][alpha.convergence==0,],2,var,na.rm=TRUE)
alpha.se.coding <- apply(results.coding[[3]][alpha.convergence==0,],2,mean,na.rm=TRUE)
alpha.coverage.coding <- apply(results.coding[[7]][alpha.convergence==0,],2,mean,na.rm=TRUE)

# pseduolikelihood
alpha.values.pl <- results.pl[[1]]
alpha.est.pl <- apply(results.pl[[1]],2,mean,na.rm=TRUE)
alpha.mcvar.pl <- apply(results.pl[[1]],2,var,na.rm=TRUE)
alpha.se.pl <- apply(results.pl[[3]],2,mean,na.rm=TRUE)
alpha.coverage.pl <- apply(results.pl[[7]],2,mean,na.rm=TRUE)

## BETA ##

# coding
beta.convergence <- results.coding[[10]]
beta.values.coding <- results.coding[[2]][beta.convergence==0,]
beta.est.coding <- apply(results.coding[[2]][beta.convergence==0,],2,mean,na.rm=TRUE)
beta.mcvar.coding <- apply(results.coding[[2]][beta.convergence==0,],2,var,na.rm=TRUE)
beta.se.coding <- apply(results.coding[[4]][beta.convergence==0,],2,mean,na.rm=TRUE)
beta.coverage.coding <- apply(results.coding[[8]][beta.convergence==0,],2,mean,na.rm=TRUE)

# pseduolikelihood
beta.values.pl <- results.pl[[2]]
beta.est.pl <- apply(results.pl[[2]],2,mean,na.rm=TRUE)
beta.mcvar.pl <- apply(results.pl[[2]],2,var,na.rm=TRUE)
beta.se.pl <- apply(results.pl[[4]],2,mean,na.rm=TRUE)
beta.coverage.pl <- apply(results.pl[[8]],2,mean,na.rm=TRUE)


############################################
############   FIGURE SETUP    #############
############################################

### MAKE DATASET FOR ALPHA ###

alphanames <- c("tau[1]","tau[2]","tau[3]",
                "rho[12]", "rho[13]", "rho[23]", 
                "nu[11]", "nu[12]", "nu[13]",
                "nu[22]", "nu[21]", "nu[23]",
                "nu[33]", "nu[31]", "nu[32]")

plotresults.alpha <- data.frame(alpha.values.coding)
colnames(plotresults.alpha) <- alphanames
plotresults.alpha.new <- cbind(melt(plotresults.alpha),1)
colnames(plotresults.alpha.new) <- c("variable","value","type")

plotresults.alpha.pl <- data.frame(alpha.values.pl)
colnames(plotresults.alpha.pl) <- alphanames
plotresults.alpha.pl.new <- cbind(melt(plotresults.alpha.pl),2)
colnames(plotresults.alpha.pl.new) <- c("variable","value","type")

plotresultsdf.alpha <- rbind(plotresults.alpha.new,plotresults.alpha.pl.new)
plottruth.alpha <- data.frame(variable = alphanames,
                        truthvalue = c(alpha.truth))

plotresultsdf.alpha %>% group_by(variable) %>%
  mutate(maxval = max(value)+.5) %>% dplyr::select(-value) %>% 
  unique() %>% as.data.frame() -> maxdf.alpha

maxdf.alpha <- data.frame(maxdf.alpha,coverage = round(c(alpha.coverage.coding,alpha.coverage.pl),2))


### MAKE DATASET FOR BETA ###
betanames <- c("beta[0]", "beta[1]", "beta[2]",
                "beta[3]", "beta[4]", "beta[5]",
                "beta[6]","beta[7]","beta[8]", "beta[9]")

plotresults.beta <- data.frame(beta.values.coding)
colnames(plotresults.beta) <- betanames
plotresults.beta.new <- cbind(melt(plotresults.beta),1)
colnames(plotresults.beta.new) <- c("variable","value","type")

plotresults.beta.pl <- data.frame(beta.values.pl)
colnames(plotresults.beta.pl) <- betanames
plotresults.beta.pl.new <- cbind(melt(plotresults.beta.pl),2)
colnames(plotresults.beta.pl.new) <- c("variable","value","type")

plotresultsdf.beta <- rbind(plotresults.beta.new,plotresults.beta.pl.new)
plottruth.beta <- data.frame(variable = betanames,
                        truthvalue = c(beta.truth))

plotresultsdf.beta %>% group_by(variable) %>%
  mutate(maxval = max(value)+.5) %>% dplyr::select(-value) %>% 
  unique() %>% as.data.frame() -> maxdf.beta

maxdf.beta <- data.frame(maxdf.beta,coverage = round(c(beta.coverage.coding,beta.coverage.pl),2))

############################################
################  FIGURES    ###############
############################################

### ALPHA ###

plotalpha <- ggplot(plotresultsdf.alpha, aes(x = factor(type), y = value, fill=factor(type))) +
  geom_boxplot() +
  pretty_plot() + labs(x = "", y = "estimates", fill = "") +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed,ncol=3) +
  geom_hline(data = plottruth.alpha, aes(yintercept = truthvalue)) +
  theme(legend.position="bottom",axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),legend.text=element_text(size=12),
        strip.text.x = element_text(size = 12)) +
  geom_text(mapping=aes(x = factor(type), y = maxval, label = coverage),
            inherit.aes = FALSE, data = maxdf.alpha, size=4,color="red") +
  scale_fill_manual(values = jdb_palette("brewer_spectra"),labels=c("Coding","PL")) 


### BETA ###

plotbeta <- ggplot(plotresultsdf.beta, aes(x = factor(type), y = value, fill=factor(type))) +
  geom_boxplot() +
  pretty_plot() + labs(x = "", y = "estimates", fill = "") +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed,ncol=3) +
  geom_hline(data = plottruth.beta, aes(yintercept = truthvalue)) +
  theme(legend.position="bottom",axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),legend.text=element_text(size=12),
        strip.text.x = element_text(size = 12)) +
  geom_text(mapping=aes(x = factor(type), y = maxval, label = coverage),
            inherit.aes = FALSE, data = maxdf.beta, size=4, color="red") +
  scale_fill_manual(values = jdb_palette("brewer_spectra"),labels=c("Coding","PL"))

##################################
############   SAVE   ############
##################################

ggsave(plotalpha, file = paste0("figures/alpha_boxplot_",as.character(samplesize),".png"), dpi = 200, width = 8, height = 10)
ggsave(plotbeta, file = paste0("figures/beta_boxplot_",as.character(samplesize),".png"), dpi = 200, width = 8, height = 10)

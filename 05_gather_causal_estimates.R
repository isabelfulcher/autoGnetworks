library(readr)
library(BuenColors)
library(tidyverse)
library(reshape2)

#args <- commandArgs(trailingOnly = TRUE)
#samplesize <- args[1]

samplesize <- "200"

#############################
#######   LOAD DATA   #######
#############################

iter <- c(1:1000)
#truth
load(paste0("input/truth_causal_estimands_",samplesize,".rda"))

#coding
file_name_coding <- paste0("output/estimates_coding_",samplesize,"_iter")

#pl
file_name_pl <- paste0("output/estimates_pl_",samplesize,"_iter")

#######################################
#######  CONSOLIDATE RESULTS    #######
#######################################

#getter function
getter <- function(iter, index) iter[[index]]

## CODING ##
readRDSx.coding <- function(xx) readRDS(paste0(file_name_coding,as.character(xx),".rds"))
allLists.coding <- lapply(iter, readRDSx.coding) #remove ones that did not converge

estimates.coding <- do.call(rbind,lapply(allLists.coding, getter, 1))[,-2]
median.boot.coding <- do.call(rbind,lapply(allLists.coding, getter, 2))[,-2]
var.coding <- do.call(rbind,lapply(allLists.coding, getter, 3))[,-2]
coverage.quantile.coding <- do.call(rbind,lapply(allLists.coding, getter, 4))[,-2]
coverage.wald.coding <- do.call(rbind,lapply(allLists.coding, getter, 5))[,-2]

coverage.coding <- apply(coverage.quantile.coding,2,mean)

## PL ##
readRDSx.pl <- function(xx) readRDS(paste0(file_name_pl,as.character(xx),".rds"))
allLists.pl <- lapply(iter, readRDSx.pl)

estimates.pl <- do.call(rbind,lapply(allLists.pl, getter, 1))[,-2]

# file list 
# x <- list.files(path = "output/", pattern = "estimates_coding_100_iter46.rds", all.files = TRUE,full.names = FALSE)
# x <- x[nchar(x)==31]
# which(c(10:99) %in% as.numeric(str_sub(str_sub(x,start=-6),end=-5)) == FALSE) + 10

############################################
############   FIGURE SETUP    #############
############################################


### MAKE DATASET FOR CODING ###
truth <- c(apply(check.truth,2,mean)[1],truth.final)[-2]
names <- c("beta(alpha)","Overall","Direct","Spillover")[-2]

plotresults <- data.frame(estimates.coding)
colnames(plotresults) <- names
plotresults.new <- cbind(melt(plotresults),1)
colnames(plotresults.new) <- c("variable","value","type")

plotresults.pl <- data.frame(estimates.pl)
colnames(plotresults.pl) <- names
plotresults.pl.new <- cbind(melt(plotresults.pl),2)
colnames(plotresults.pl.new) <- c("variable","value","type")

plotresultsdf <- rbind(plotresults.new,plotresults.pl.new)
plottruth <- data.frame(variable = names,
                        truthvalue = c(truth))

plotresultsdf %>% group_by(variable) %>%
  mutate(maxval = max(value)+.05) %>% dplyr::select(-value) %>% 
  unique() %>% as.data.frame() -> maxdf

maxdf <- data.frame(maxdf,coverage = round(c(coverage.coding,rep(NA,length(coverage.coding))),2))

############################################
################  FIGURES    ###############
############################################

plot <- ggplot(plotresultsdf, aes(x = factor(type), y = value, fill=factor(type))) +
  geom_boxplot() +
  pretty_plot() + labs(x = "", y = "estimates", fill = "") +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed,ncol=3) +
  geom_hline(data = plottruth, aes(yintercept = truthvalue)) +
  theme(legend.position="bottom",axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),legend.text=element_text(size=12),
        strip.text.x = element_text(size = 12)) +
  geom_text(mapping=aes(x = factor(type), y = maxval, label = coverage),
            inherit.aes = FALSE, data = maxdf, size=4, color="red") +
  scale_fill_manual(values = jdb_palette("brewer_spectra"),labels=c("Coding","PL"))

############################################
################  TABLES  ##################
############################################

prop_bias <- function(x) { (x - truth)/truth }



table.coding <- round(data.frame(truth=truth,
                           mean_est=apply(estimates.coding,2,mean), 
                         #  prop_bias = apply(apply(estimates.coding,1,prop_bias),1,mean),
                           bias = apply(estimates.coding,2,mean)-truth,
                           variance = apply(var.coding,2,mean),
                           mc_variance = apply(estimates.coding,2,var), 
                           coverage = coverage.coding),3)


table.pl <-  data.frame(truth=round(truth,3),
                       mean_est= round(apply(estimates.pl,2,mean),3), 
                      # prop_bias = round(apply(apply(estimates.pl,1,prop_bias),1,mean),3),
                       bias = round(apply(estimates.coding,2,mean)-truth,3),
                       variance = NA,
                       mc_variance = round(apply(estimates.pl,2,var),3), 
                       coverage = NA)

samplesize
table.coding
table.pl


##################################
############   SAVE   ############
##################################

file_name_output <- paste0("figures/causal_estimates_coding_table_",samplesize,".txt")
write.table(table.coding,file_name_output, quote = FALSE, sep = "\t")

file_name_output <- paste0("figures/causal_estimates_pl_table_",samplesize,".txt")
write.table(table.pl,file_name_output, quote = FALSE, sep = "\t")

ggsave(plot, file = paste0("figures/causal_estimates_boxplot_",as.character(samplesize),".png"), dpi = 200, width = 10, height = 6)

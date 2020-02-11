#Read in Data
Absolute_supp <- read.csv('~/Desktop/Research/CCLE/AbsoluteSupplement1.csv')

#Create new dataframe containing only the data from cell lines
Absolute_supp_CL <- Absolute_supp[Absolute_supp$TUMOR_TYPE == 'Cell-line',]

#Delete the hyphens in the cell line names to enable merging with expression data
Absolute_supp_CL$SAMPLE_NAME <- gsub('-', '', Absolute_supp_CL$SAMPLE_NAME)

#Write out DataSet to edit in BBedit
write.csv(Absolute_supp_CL, '~/Desktop/Research/CCLE/SortedAbsoluteSupp.csv')

#Read in DataSet edited in BBedit
Absolute_supp_CL2 <- read.csv('~/Desktop/Research/CCLE/SortedAbsoluteSupp2.csv')

#Remove unwanted columns
Absolute_supp_CL2$ARRAY <- NULL
Absolute_supp_CL2$MATCHED_PLATE_NORMAL_ARRAY <- NULL
Absolute_supp_CL2$TUMOR_NORMAL <- NULL
Absolute_supp_CL2$TISSUE_SITE <- NULL
Absolute_supp_CL2$OSH_DIAGNOSIS <- NULL
Absolute_supp_CL2$LINEAGE <- NULL
Absolute_supp_CL2$COHORT <- NULL
Absolute_supp_CL2$PLATE <- NULL
Absolute_supp_CL2$PLATFORM <- NULL
Absolute_supp_CL2$X <- NULL

#rename cell line variable label to allow merging with expression data
names(Absolute_supp_CL2)[1] <- 'line'

#Merge Edited Absolute supplement and Expression datasets
#See 'Expression_Absolute_HighLow.R' for details on expression dataset
Absolute_supp_Expression <- merge(Absolute_supp_CL2, Expression_Methylation_3, by = 'line')

#Make the a DataFrame with numeric values only and 'line' as title
Absolute_supp_Expression_2 <- Absolute_supp_Expression[,-(1:6)]
rownames(Absolute_supp_Expression_2) <- (Absolute_supp_Expression[,1])

Absolute_supp_Expression_4 <- Absolute_supp_Expression[,-(1:4)]
Absolute_supp_Expression_4 <- Absolute_supp_Expression_4[,-(2)]
rownames(Absolute_supp_Expression_4) <- (Absolute_supp_Expression[,1])

#Make the DataFrame without NA values
Absolute_supp_Expression_3 <- na.omit(Absolute_supp_Expression_2)

Absolute_supp_Expression_4 <- na.omit(Absolute_supp_Expression_4)


#Function to make a new dataframe with mean, sd, se, and sample size
#x is dataframe and variable of interest
#y is dataframe of interest
#a is low value, ex: .25 for lowest 25%
#b is high value
#c is variable to measure
SplitExpression3 <- function(x,y,a,b){
  q1 <- unname(quantile(x, a))
  q2 <- unname(quantile(x, b))
  Low <- y[x <= q1,]
  High <- y[x >= q2,]
  
  AvgLow_scgf <- mean(Low$Subclonal.genome.fraction)
  AvgHigh_scgf <- mean(High$Subclonal.genome.fraction)
  sdLow_scgf <- sd(Low$Subclonal.genome.fraction)
  sdHigh_scgf <- sd(High$Subclonal.genome.fraction)
  #seLow_scgf <- sdLow_scgf/sqrt(nrow(Low$Subclonal.genome.fraction))
  #seHigh_scgf <- sdHigh_scgf/sqrt(nrow(High$Subclonal.genome.fraction))
  seLow_scgf <- sdLow_scgf/sqrt(nrow(Low[,1, drop=FALSE]))
  seHigh_scgf <- sdHigh_scgf/sqrt(nrow(High[,1, drop=FALSE]))
  nLow <- nrow(Low)
  nHigh <- nrow(High)
  
  variable <- c('low', 'high')
  value <- c(AvgLow_scgf, AvgHigh_scgf)
  sd <- c(sdLow_scgf, sdHigh_scgf)
  se <- c(seLow_scgf, seHigh_scgf)
  n <- c(nLow, nHigh)
  
  Split2 <- data.frame(variable, value, sd, se, n)
}

library(ggplot2)
theme_set(theme_classic())

KDM4A_supp_3 <- SplitExpression3(Absolute_supp_Expression_3$KDM4A, Absolute_supp_Expression_3, .25, .75)

ggplot(KDM4A_supp_3, aes(variable, value)) + 
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') + 
  geom_errorbar(aes(ymin=value - se, ymax=value + se), width = 0.25, position = position_dodge(.9))

ggsave(filename = './KDM4A_25.tiff', scale = 1, width = 7.5, height = 7.5, units = 'cm', dpi = 300)


PHF8_supp_3 <- SplitExpression3(Absolute_supp_Expression_3$PHF8, Absolute_supp_Expression_3, .25, .75)

ggplot(PHF8_supp_3, aes(variable, value)) + 
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') + 
  geom_errorbar(aes(ymin=value - se, ymax=value + se), width = 0.25, position = position_dodge(.9))
  
ggsave(filename = './PHF8_25.tiff', scale = 1, width = 7.5, height = 7.5, units = 'cm', dpi = 300)



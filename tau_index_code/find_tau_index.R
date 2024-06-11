#read in median expression data 
median_expression <- read.csv("median_expression.csv")

#create min-max normalisation function - the tau paper detailed that expression 
#values were normalised to the maximum value for a given tissue per gene. Min-max
#normalisation seems appropriate here, as this normalises to the largest expression
#value and ensures expression values are on a scale from 0 to 1 
min_max_normalisation <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}

#use apply function to convert data frame into min-max normalised data
normalised_expression <- as.data.frame(t(apply(median_expression, 
                                               1, min_max_normalisation)))

#create tau index function 
tau_index <- function(x){
  return(sum(1-x)/length(x))
}

#Add Tau index column to data
normalised_expression$Tau_index <- apply(normalised_expression, 1, tau_index)


                                       
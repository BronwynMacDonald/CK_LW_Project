# Bronwyn MacDonald
# May 9 2017

# Fitting von Bertalanffy Model to Chinook data
# For evaluation of size limits in ocean fisheries 

# Things to NOTE:
# if CWT==TRUE it estiamtes the age composition using the age LVB model
# if CWT==FALSE it uses the non-age structured model
# If POOLED == TRUE you do not mneed to hcange the assignments file, the code sets the POP_ID and REF_ID of MIddle Shuswap equal to that of 
# Lower Shuswap and no further adjustments have to be made
rm(list=ls())

Pool_Shu=TRUE   # Pool the Lower and Middle Shuswap data 



# =============== Step 1: Background =====================
# Code that installs packages if need be - taken from Forecast.R
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

# Install packages
usePackage("lattice")
usePackage("superdiag")
usePackage("FSA")
usePackage("dplyr")
usePackage("reshape")
usePackage("ggplot2")
usePackage("mcmcplots")
usePackage("nlstools")
usePackage("RColorBrewer")
usePackage("abind")

# Source functions and models
source("LVB Models.R")

# If data have not been formatted and transformed must run the clean data code
# source("Data Cleanup_2.R)

# Read in clean data
AL_DATA <- read.csv("CK Age-Length Data.csv")
# Read in reference stocks for age correction 
assignments <- read.csv("theta matrices.csv")

# =================== Step 2: Correct scale age bias

# Get paired CWT-Scale age data to estimate true age probabilities based on scale ages
# When running the age error models must define the Age and the True Age
  Both_Data <- Get_Both_Data(AL_DATA) 
  if(Pool_Shu == TRUE)  {
    Both_Data <- Get_Both_Pooled(AL_DATA)
  }
  Both_Plot <- as.data.frame(Both_Data)


# Fit the CWT-Scale Age relationship (and lVB model) to data that is paired
  CWT_Only <- filter(Both_Plot, !is.na(True))

# Fit the LVB age proportion model to get thetas
  CWT_Model <- Run_Model(Data=CWT_Only,Bayes_Data=CWT_Only, LOG=TRUE,POP=TRUE ,nboots=300, Freq=TRUE, CWT=TRUE, thetas=FALSE)
  write.csv(CWT_Model$thetas, "thetas.csv")

# If adjusted Nicola is being used for Spring and Summer 52 then must create thetas 
  if(nrow(filter(assignments, Ref_CWT =="Nicola_Adj")) >0) thetas.out <- Nicola_Calc(CWT_Only, CWT_Model$thetas) else(thetas.out <- thetas)

#Lookup table of reference IDs for each CWT stock to be applied to scale aged populations
  Lookup <- data.frame(Pops = unique(CWT_Only$Pop), ID=CWT_Only$POP_ID[match(unique(CWT_Only$Pop), CWT_Only$Pop)]) # if POOLED both Shu should = 2
  if(nrow(filter(assignments, Ref_CWT =="Nicola_Adj")) >0){
    levels(Lookup$Pops) <- c(levels(Lookup$Pops), "Nicola_Adj")
    Lookup <- rbind(Lookup, c("Nicola_Adj", dim(thetas.out)[3]))
  }
  assignments$REF_ID <- Lookup$ID[match(assignments$Ref_CWT,Lookup$Pops)]
  
# Pull Scale aged only data
  Scale_Data <- Get_Scale_Data(AL_DATA) 
  Scale_plot <- as.data.frame(Scale_Data)

# Assign reference IDs to scale aged populations  
  Scale_plot$REF_ID <- as.numeric(assignments$REF_ID[match(Scale_plot$Pop,assignments$Stock)]) 
  
# Estimate true ages
  Est_True <- AgePred(Data=Scale_plot, alpha=thetas.out) # Runs the age calculation using pre-determined age matrix from CWT fit
                                                         # thetas are an array with nPop slices of 6X6 tables of probabilities (median probabilities)
                                                         # Must input data used to create theta matrix
  #Est_True_2 <- AgePred_noprior(Data=Scale_plot, theta=thetas.out)
# Age matrix has true age estimates as well as theta estimates
  Est.age.mat <- data.frame(as.matrix(Est_True))
  Pred.Age <-  dplyr::select(Est.age.mat, starts_with("True"))
  pred.thetas <-  dplyr::select(Est.age.mat, starts_with("theta"))
  #True Age (age predictions are normal)
  pred.true.median <- apply(Pred.Age, 2, FUN=quantile, probs=0.5)
  #thetas
  pred.thetas.median <- structure(.Data = apply(pred.thetas, 2,  FUN=quantile, probs=0.5), 
                              .Dim=c(6,6,dim(thetas.out)[3]))
  pred.thetas.mode <- structure(.Data = apply(pred.thetas, 2,  FUN=Mode), 
                                .Dim=c(6,6,dim(thetas.out)[3]))

# Input the estimated true ages into the data file!
  Scale_plot$True <- pred.true.median
  print(table(Scale_plot$True, Scale_plot$Age, Scale_plot$Pop))

# ====================== Step 3: Select Data Aggregation and Run  the LVB with the estimated true ages ==========

 table(Scale_plot$Pop)
 table(Scale_plot$MU)
# Run as POPULATIONS
# If running populations must remove those with too little data (<10 data points)
 Pop_Data <- filter(Scale_plot, Pop %in%(levels(Scale_plot$Pop)[tabulate(Scale_plot$Pop)>10]) )
 # Re-number POP_ID so is consecutive now that some have been excluded
 Pop_Data$POP_ID=match(Pop_Data$Pop, unique(Pop_Data$Pop))
 Pop_Model <- Run_Model(Data=Pop_Data,Bayes_Data=Pop_Data, LOG=TRUE,POP=TRUE ,nboots=300, Freq=FALSE, CWT=FALSE, TrueEst=TRUE) 
 
 # Plot Pars
 Boxplot(Pop_Model$Model, Pop_Data, "Pop_Model")
 Single_Plot(Pop_Model$Model, Pop_Data, "Pop_Model")
# Run as MU

# Run with Estimated true age
  Est_Model <- Run_Model(Data=Scale_plot,Bayes_Data=Scale_plot, LOG=TRUE,POP=TRUE ,nboots=300, Freq=FALSE, CWT=FALSE, TrueEst=TRUE) 


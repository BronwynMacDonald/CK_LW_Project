
rm(list=ls())

# Data Clean-up
library(plyr)
library(dplyr)
library(ggplot2)
library(R2jags)


Raw_data <- read.csv("CN master flat file June 9.csv")  # 284,249 samples 


Age_data <- filter(Raw_data, Standard.G.R.scale.ages != "" | Standard.G.R.scale.ages != "2M" | Standard.G.R.scale.ages != "3M" | CWT.age != "") # 56,483 samples 
Age_length <- filter(Age_data, FORK.Length != "" | POH != "") # 56,214


# Pull out date information
  LIVE_dates <- as.Date(as.character(Age_length$Live.Date..YYYY.MM.DD), format = "%d/%m/%Y")
  Age_length$LIVE_Month <- format(LIVE_dates, "%m")
  Age_length$LIVE_Day <- format(LIVE_dates, "%d")
  REC_dates <- as.Date(as.character(Age_length$REC.Date....YYYY.MM.DD), format = "%d/%m/%Y")
  Age_length$REC_Month <- format(REC_dates, "%m")
  Age_length$REC_Day <- format(REC_dates, "%d")

  
  
# Convert numbers to numeric and extract Ocean and Stream-type POH
Age_length$POH <- as.numeric(as.character(Age_length$POH))*10 
Age_length$FORK.Length <- as.numeric(as.character(Age_length$FORK.Length))

# Create dataframe
AL_DATA <- select(Age_length, Year, MU=MU.Name, Type=Stream.or.Ocean.Type, CU=CU, CTC=CTC.Model.Stock.Name, Pop=Population, Stream=Stream.Name, LIVE_Month, LIVE_Day, LIVE_Sex=Live..Sex, 
                  LIVE_FL = FORK.Length, LIVE_Adipose=Live.Adipose, REC_Month, REC_Day, REC_Sex=REC.Sex, REC_Adipose=Rec.Adipose, Recovered=recovered.y.n.,  
                  Sampled=Sampled..Y.N., POH, CWT_Age=CWT.age, CWT_BY=CWT.brood.year, Scale_Age=Standard.G.R.scale.ages
                  )
rm(list=c("Raw_data", "Age_length"))
# Remove data that does not have POH Length
AL_DATA <- AL_DATA %>%  filter(!is.na(POH)) # 55,875
AL_DATA$Est_FL <- AL_DATA$Est_FL_lower <- AL_DATA$Est_FL_upper <-NA 
 

Ocean.POH <-  AL_DATA %>% filter(Type == "Ocean")  %>%           # 43,608
                          select(POH)  


Stream.POH  <-  AL_DATA %>% filter(Type == "Stream") %>%         # 12,267  
                            select(POH)    


# Run the POH to Fork Length conversion
source("POH_FL Conversion.R")


#  =========== Calculate the Ocean-type Fork Lengths ==================== #

regression_Ocean <- Regression(POH_FL_data_2008_4_1, new_data =Ocean.POH[,1])
#regression_Ocean$BUGSoutput$summary
Ocean_Lower_2.5 <- select(data.frame(as.list(regression_Ocean$BUGSoutput$summary[,"2.5%"])), -alpha, -beta, -deviance, -tau)
Ocean_50p <- select(data.frame(as.list(regression_Ocean$BUGSoutput$summary[,"50%"])), -alpha, -beta, -deviance, -tau)
Ocean_Upper_97.5 <- select(data.frame(as.list(regression_Ocean$BUGSoutput$summary[,"97.5%"])), -alpha, -beta, -deviance, -tau)


AL_DATA$Est_FL[AL_DATA$Type == "Ocean"] <- as.vector(unlist(Ocean_50p))
AL_DATA$Est_FL_lower[AL_DATA$Type == "Ocean"] <- as.vector(unlist(Ocean_Lower_2.5))
AL_DATA$Est_FL_upper[AL_DATA$Type == "Ocean"] <- as.vector(unlist(Ocean_Upper_97.5))


#  =========== Calculate the Stream-type Fork Lengths ==================== #

regression_Stream <- Regression(POH_FL_data_1981_5_2, new_data = Stream.POH[,1])
#regression_Stream$BUGSoutput$summary
Stream_Lower_2.5 <- select(data.frame(as.list(regression_Stream$BUGSoutput$summary[,"2.5%"])), -alpha, -beta, -deviance, -tau)
Stream_50p <- select(data.frame(as.list(regression_Stream$BUGSoutput$summary[,"50%"])), -alpha, -beta, -deviance, -tau)
Stream_Upper_97.5 <- select(data.frame(as.list(regression_Stream$BUGSoutput$summary[,"97.5%"])), -alpha, -beta, -deviance, -tau)

AL_DATA$Est_FL[AL_DATA$Type == "Stream"] <- as.vector(unlist(Stream_50p))
AL_DATA$Est_FL_lower[AL_DATA$Type == "Stream"] <- as.vector(unlist(Stream_Lower_2.5))
AL_DATA$Est_FL_upper[AL_DATA$Type == "Stream"] <- as.vector(unlist(Stream_Upper_97.5))


# Look at the predicted values 
qplot(POH, Est_FL, data=AL_DATA, colour=Type)

output_Ocean <- as.mcmc(regression_Ocean)
output_Stream <- as.mcmc(regression_Stream)



# Compare alphas
alpha_pos <- which(varnames(output_Ocean)=="alpha")
boxplot(unlist(output_Stream[,alpha_pos]), unlist(output_Ocean[,alpha_pos]), main="Alphas", names=c("Stream-type","Ocean-type"))

# Compare betas
beta_pos <- which(varnames(output_Ocean)=="beta")
boxplot(unlist(output_Stream[,beta_pos]), unlist(output_Ocean[,beta_pos]), main="Betas", names=c("Stream-type","Ocean-type"))



# ======================== Look at age data ================================ #

# CWT Data = 13,650 of 55,875
# Scale age data = 6726 empty

AL_DATA$Scale_Agg_Age <- NA  
AL_DATA$Scale_Agg_Age[AL_DATA$Scale_Age=="21" ] <- 2   # 4,593
AL_DATA$Scale_Agg_Age[AL_DATA$Scale_Age=="31" | AL_DATA$Scale_Age=="32" ] <- 3   # 10,842
AL_DATA$Scale_Agg_Age[AL_DATA$Scale_Age=="41" | AL_DATA$Scale_Age=="42"] <- 4  #29,215
AL_DATA$Scale_Agg_Age[AL_DATA$Scale_Age=="51" | AL_DATA$Scale_Age=="52"] <- 5  # 4,388
AL_DATA$Scale_Agg_Age[AL_DATA$Scale_Age=="61" | AL_DATA$Scale_Age=="62"] <- 6  # 110

# Paired_data <- AL_DATA %>%  select(CWT_Age, Scale_Agg_Age) %>%
#                             filter(CWT_Age != "" & !is.na(Scale_Agg_Age))   # 6,924 samples
# 
# Paired_data <- data.frame(lapply(Paired_data, as.numeric))
# Paired.t <- table(Paired_data)
# N.Samples <- as.vector(Paired.t)
# N.Samples[N.Samples==0]<-NA
# qplot( rep(seq(2,5),5),  rep(2:6, each=4), cex=N.Samples, colour=N.Samples, ylim=c(2,5), xlab="CWT Age", ylab="Scale Age")
# print(Paired.t)
# 
# # error in scale ages is 0.0702, 7% of samples; 2.3% over-estimates of age; 4.65% under-estimates
# # 
# Paired_stocks <-  AL_DATA %>%   filter(CWT_Age != "" & !is.na(Scale_Agg_Age))   # 6,924 samples
# Harrison <- Paired_stocks %>%   filter(Pop=="Harrison")
# LowerShu <- Paired_stocks %>%   filter(Pop=="Shuswap R (Lower)")                                
# Nicola <- Paired_stocks %>%   filter(Pop=="Nicola")
# MiddleShu <- Paired_stocks %>%   filter(Pop=="Shuswap R (Middle)")
# Chilliwack <- Paired_stocks %>%   filter(Pop=="Chilliwack")

#AL_DATA_DF <- data.frame(lapply(AL_DATA, as.character), stringsAsFactors=FALSE)
write.csv(AL_DATA, "CK Age-Length Data.csv")


# pred<-data.frame( CWT_Age=Harrison$CWT_Age, predict(lm(Scale_Agg_Age~CWT_Age, data=Harrison), interval="confidence") )
# 
# 
# 
# require(ggplot2)
# ggplot(pred, aes(x = CWT_Age, y = fit)) +
#   geom_point(aes(x=CWT_Age, y=fit), size = 2,colour="dodgerblue") +
#   geom_errorbar(aes(ymax = upr, ymin = lwr),size=0.5, width=0.15,color="dodgerblue") +
#   geom_abline(slope=1, intercept=0, color="dodgerblue3") +
#   ylab("Fit Scale Age") +
#   xlab("CWT Age") +
#   scale_y_continuous(minor_breaks = NULL) +
#   scale_x_continuous(minor_breaks = NULL)


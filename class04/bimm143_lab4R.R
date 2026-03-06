setwd("~/Desktop/bimm 143")

source("http://thegrantlab.org/misc/cdc.R")
View(cdc)
head(cdc$height)

#question 1
tail(cdc$weight, 20)

#libraries
library(ggplot2)
library(vegan)

#question 2
#dry mix (weight depends on height)
plot<-ggplot(cdc, mapping = aes(height, weight)) + geom_point(col= "darkgreen") + theme_bw() + labs(title = "Weight vs Height Plot", cex =4)

#question 5 
#do they appear to be correlated? yes

#question 4
cor(cdc$height, cdc$weight)

#question 5
height_m <- cdc$height * 0.0254
weight_kg <- cdc$weight * .454
weight_m

#question 6
BMI <- (weight_kg)/(height_m^2)
plot(cdc$height, BMI)

#question 7
cor(cdc$height, BMI)

#question 8 
sum(BMI >= 30)

#question 9
height_data <- cdc[1:100, 5]
height_data
weight_data <- cdc[1:100, 6]
weight_data
plot(height_data, weight_data)

#question 10
gender_data <- cdc[which(cdc$gender=="m" & BMI >30),]
gender_data

gender_data_2 <- subset(gender_data, gender = "m" | BMI >= 30)


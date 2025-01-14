# install.packages("data.table")
# install.packages("raster")
# install.packages("randomForest")
# install.packages("lattice")
# install.packages("RColorBrewer")
# install.packages("PresenceAbsence")
library(raster)

avi_dat <- read.table('data/Data_SwissBreedingBirds.csv', header=T, sep=',')

summary(avi_dat)

avi_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'std', 'rad', 'blockCV_tile')

avi_df <- data.frame(avi_dat)[avi_cols]

summary(avi_df)

# Please note that you have to set download=T if you haven't downloaded the data before:
# NO LONGER WORKS!
# bio_curr <- getData('worldclim', var='bio', res=0.5, lon=5.5, lat=45.5, path='data')[[c(2,5,14)]]

# Mark / Hanna's code to get data 
# requires "terra" package 
# bio_curr <- worldclim_country("Switzerland", var="bio", path="data", download=T)

# Please note that you have to set download=T if you haven't downloaded the data before:
# didn't get to the future data 
# bio_fut <- getData('CMIP5', var='bio', res=0.5, lon=5.5, lat=45.5, rcp=45, model='NO', year=50, path='data', download=F)[[c(2,5,14)]]

bio_curr_file <- "data/CHE_wc2.1_30s_bio.tif"
bio_curr <- stack(bio_curr_file)
plot(bio_curr)

bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
ch_ext <- c(5, 11, 45, 48)
bio_curr <- crop(bio_curr, ch_ext)

# Re-project to Swiss coordinate system and clip to Swiss political boundary
bio_curr <- projectRaster(bio_curr, bg)
bio_curr <- resample(bio_curr, bg)
bio_curr <- mask(bio_curr, bg)

names(bio_curr) # but we want c('bio_2', 'bio_5', 'bio_14')
#subset to variables of interest 
bio_curr <- bio_curr[[c("wc2.1_30s_bio_2","wc2.1_30s_bio_5","wc2.1_30s_bio_14")]]

plot(bio_curr)


m_glm <- glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=avi_df)
summary(m_glm)

library(mecofun)

# Names of our variables:
pred <- c('bio_2', 'bio_5', 'bio_14')

# We want three panels next to each other:
par(mfrow=c(1,3)) 

# Plot the partial responses
partial_response(m_glm, predictors = avi_df[,pred])



library(RColorBrewer)
library(lattice)

# We prepare the response surface by making a dummy data set where two predictor variables range from their minimum to maximum value, and the remaining predictor is kept constant at its mean:
xyz <- data.frame(expand.grid(seq(min(avi_df[,pred[1]]),max(avi_df[,pred[1]]),length=50), seq(min(avi_df[,pred[2]]),max(avi_df[,pred[2]]),length=50)), mean(avi_df[,pred[3]]))
names(xyz) <- pred

# Make predictions
xyz$z <- predict(m_glm, xyz, type='response')
summary(xyz)

# Make a colour scale
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# plot 3D-surface
wireframe(z ~ bio_2 + bio_5, data = xyz, zlab = list("Occurrence prob.", rot=90),
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          main='GLM', xlab='bio_2', ylab='bio_5', screen=list(z = 120, x = -70, y = 3))

# Plot inflated response curves:
par(mfrow=c(1,3)) 
inflated_response(m_glm, predictors = avi_df[,pred], method = "stat6", lwd = 3, main='GLM') 

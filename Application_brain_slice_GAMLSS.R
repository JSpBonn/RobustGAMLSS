# Application brain data:

library(gamair)
data(brain) # 1567 observations (no exclusions)

#setwd("//folderapplication")
source(Robust_gamboostLSS_Families.R) # loading some packages and the robust Gamma families object for GAMLSS gradient boosting


stopping=2000

# MAD non cyclic
set.seed(123) # ordinary gamma distribution for GAMLSS with boosting
gam3 <- gamboostLSS(formula=medFPQ~bspatial(Y,X), data = brain[,1:3],method = "noncyclic",control = boost_control(mstop=stopping),families=GammaLSS(stabilization = "MAD"))


c_value0_05 <- c_generate_Gamma(brain$medFPQ)
set.seed(123) # robust gamma distribution for GAMLSS with boosting
gam4 <-gamboostLSS(formula = medFPQ~bspatial(Y,X),data= brain[,1:3], method = "noncyclic", control = boost_control(mstop=stopping), families = robust_GammaLSS(stabilization = "MAD", rob=c_value0_05)) 

set.seed(1234)
cvr_gam3<- cvrisk(gam3)
m_gam3 <- mstop(cvr_gam3) # optimal mstop = 1480
# m_gam3 = 1480 # shortcut 

set.seed(1234)
cvr_gam4<- cvrisk(gam4)
m_gam4 <- mstop(cvr_gam4) # optimal mstop = 197 
# m_gam4 = 197 # shortcut 

##### parametrisation of gamboostLSS package for gamma distribution 
##### additive predictors and early stopping (predictions based on early stopping)

pred_mu <- log(predict(gam3[m_gam3],newdata = brain[,1:3],parameter = "mu",type="response")) # eta_theta_mu, non robust
pred_sig <- log(predict(gam3[m_gam3],newdata =brain[,1:3],parameter = "sigma",type="response"))  # eta_theta_sigma, non robust
pred_mu2 <- log(predict(gam4[m_gam4],newdata = brain[,1:3],parameter = "mu",type="response")) # eta_theta_mu, robust
pred_sig2 <- log(predict(gam4[m_gam4],newdata = brain[,1:3],parameter = "sigma",type="response")) # eta_theta_sigma, robust
# same as leaving type = "response" away not not taking the logarithm, but for different parametrisations it might be useful (see below)


##### different parametrisation of a gamma distribution (expection value = mu, variance = mu^2*sigma^2, as in Aeberhard et al., 2021, doi: 10.1007/s11222-020-09979-x)
##### additive predictors and predictions on converged coefficients

#pred_mu <- log(predict(gam3,newdata = brain[,1:3],parameter = "mu",type="response"))# eta_theta_mu, non robust
#pred_sig <-log(1/predict(gam3,newdata = brain[,1:3],parameter = "sigma",type="response")^0.5) # eta_theta_sigma, non robust
#pred_mu2 <- log(predict(gam4,newdata = brain[,1:3],type="response",parameter = "mu")) # eta_theta_mu, robust
#pred_sig2 <- log(1/predict(gam4,newdata = brain[,1:3],type="response",parameter = "sigma")^0.5) # eta_theta_sigma, robust



library(gridExtra)
library(lattice)

# non robust
data_example <-  data.frame(brain[,1:2],pred_mu)

levelplot(pred_mu ~ Y * X, data_example,main="Eta_mu non robust",xlab=expression(paste("Voxel horizontal axis (x"[2],")",sub="")) ,ylab=expression(paste("Voxel horizontal axis (x"[1],")",sub="")),
          panel = panel.levelplot,cuts = 80)


data_example2 <-  data.frame(brain[,1:2],pred_sig)

levelplot(pred_sig ~ Y * X, data_example2,main="Eta_sigma non robust",xlab=expression(paste("Voxel horizontal axis (x"[2],")",sub="")) ,ylab=expression(paste("Voxel horizontal axis (x"[1],")",sub="")),
          panel = panel.levelplot,cuts = 80)


# robust
data_example3 <-  data.frame(brain[,1:2],pred_mu2)

levelplot(pred_mu2 ~ Y * X, data_example3,main="Eta_mu robust",xlab=expression(paste("Voxel horizontal axis (x"[2],")",sub="")) ,ylab=expression(paste("Voxel horizontal axis (x"[1],")",sub="")),
          panel = panel.levelplot,cuts = 80)


data_example4 <-  data.frame(brain[,1:2],pred_sig2)

levelplot(pred_sig2 ~ Y * X, data_example4,main="Eta_sigma robust",xlab=expression(paste("Voxel horizontal axis (x"[2],")",sub="")) ,ylab=expression(paste("Voxel horizontal axis (x"[1],")",sub="")),
          panel = panel.levelplot,cuts = 80)


# take care: 
# making them comparable within a parameter/additive predictor is additional work in R
# produce a grid, where everything is NA except the brain region (you can make it smoother with a tighter grid/lattice and predict these values from the estimated model, too)
# you can also smooth the boundaries of the brain slice to generate a prettier plot (cf. in Speller et al.) 










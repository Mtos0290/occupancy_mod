###################################################################################
## Analisis de ocupacion A. ordinarium 60 sitios temporada de secas             ###
###################################################################################

R

library(unmarked)
library(plotrix)
library(AICcmodavg)
library(usdm)
library(corrplot)




###### Importacion y manipulacion de la base de datos###############

amby <- read.table("ambystoma_datos_b.csv", sep= ",", header = T)

summary(amby)

names(amby)

############estimacion de los promedios de variables#####

amby$trans <- apply(amby[,45:46], 1, sum) ##proporción de uso de suelo correspondiente a cultivos y pastizal inducido

summary(amby$trans)

amby$mean_temp <- apply(amby[,8:10], 1, mean)

amby$mean_ph <-  apply(amby[,14:16], 1, mean)

amby$mean_cond <-  apply(amby[,20:22], 1, mean)

amby$mean_OD <-  apply(amby[,26:28], 1, mean)

###historiales de detección####

amby.y <- amby[,2:4]


#####transfomacion a logaritmo de coductividad#####

amby$lcond <- log(amby$mean_cond)

names(amby)

###########estandarización de covariables#######

amby.cov <- data.frame(scale(amby[c(38,47:50,52,53)]))

names(amby.cov)

##########verificamos el grado de correlación de las variables consideradas####

vif(amby.cov)

corrplot(cor(amby.cov), method="number", type = "upper")


#####observaciones para los modelos ####
amby.obs<-list(temp=temp, prof=prof, cobr=cobr)

####base de datos para modelos####
ambyUMF <- unmarkedFrameOccu(	y = amby.y, siteCovs = amby.cov)

summary(ambyUMF)

####gráficos de detecciones#####

plot(ambyUMF, col.regions = palette(gray(2:1/3)))


############modelos de ocupación temporada de secas##################

########################vegetación nativa######################

psi.nat500_p. <- occu(~1 ~nat_500, data = ambyUMF)

psi.nat500_temp_p. <- occu(~1 ~nat_500 + mean_temp, data = ambyUMF)

psi.nat500_lcond_p. <- occu(~1 ~nat_500 + lcond, data = ambyUMF)

psi.nat500_pH_p. <- occu(~1 ~nat_500 + mean_ph, data = ambyUMF)

psi.nat500_OD_p. <- occu(~1 ~nat_500 + mean_OD, data = ambyUMF)

psi.nat500_elev_p. <- occu(~1 ~nat_500 + elev, data = ambyUMF)

########################transformado##########################
psi.trans500_p. <- occu(~1 ~trans, data = ambyUMF)

psi.trans500_temp_p. <- occu(~1 ~trans + mean_temp, data = ambyUMF)

psi.trans500_lcond_p. <- occu(~1 ~trans + lcond, data = ambyUMF)

psi.trans500_pH_p. <- occu(~1 ~trans + mean_ph, data = ambyUMF)

psi.trans500_OD_p. <- occu(~1 ~trans + mean_OD, data = ambyUMF)

psi.trans500_elev_p. <- occu(~1 ~trans + elev, data = ambyUMF)

#########################variables locales########################

psi.temp_p. <- occu(~1 ~ mean_temp, data=ambyUMF)

psi.lcond_p. <- occu(~1 ~ lcond, data = ambyUMF)

psi.pH_p. <- occu(~1 ~ mean_ph, data = ambyUMF)

psi.OD_p. <- occu(~1 ~ mean_OD, data = ambyUMF)

psi.elev_p. <- occu(~1 ~ elev, data = ambyUMF)

psi.elev_OD_p. <- occu(~1 ~ elev + mean_OD, data = ambyUMF)

psi.temp_OD_p. <- occu(~1 ~ mean_temp + mean_OD, data = ambyUMF)

psi.lcond_OD_p. <- occu(~1 ~ lcond + mean_OD, data = ambyUMF)

nulo <- occu(~1 ~1, data = ambyUMF)

psi.elev_2_p. <- occu(~1 ~ I(elev^2), data = ambyUMF)

psi.pH_lcond_p. <- occu(~1 ~mean_ph + lcond, data = ambyUMF)

#################selección modelos #######
##############lista de modelos############


Cand.models <- list(psi.nat500_p., psi.nat500_temp_p., 
                    psi.nat500_lcond_p., psi.nat500_pH_p.,
                    psi.nat500_OD_p., psi.nat500_elev_p.,
                    psi.trans500_p.,
                    psi.trans500_temp_p.,
                    psi.trans500_lcond_p.,
                    psi.trans500_pH_p.,
                    psi.trans500_OD_p.,
                    psi.trans500_elev_p.,
                    psi.temp_p., psi.elev_OD_p.,
                    psi.lcond_p., 
                    psi.pH_p., 
                    psi.OD_p., 
                    psi.elev_p., 
                    psi.temp_OD_p.,
                    psi.lcond_OD_p.,
                    nulo, psi.elev_2_p.,
                    psi.pH_lcond_p.)

##############nombres de candidatos###########

mod.names<- c("psi.nat500_p.", "psi.nat500_temp_p.", 
              "psi.nat500_lcond_p.", "psi.nat500_pH_p.",
              "psi.nat500_OD_p.", "psi.nat500_elev_p.", 
              "psi.trans500_p.",
              "psi.trans500_temp_p.",
              "psi.trans500_lcond_p.",
              "psi.trans500_pH_p.",
              "psi.trans500_OD_p.",
              "psi.trans500_elev_p.",
              "psi.temp_p.", "psi.elev_OD_p.",
              "psi.lcond_p.", 
              "psi.pH_p.", 
              "psi.OD_p.", 
              "psi.elev_p.", 
              "psi.temp_OD_p.",
              "psi.lcond_OD_p.",
              "nulo", "psi.elev_2_p.",
              "psi.pH_lcond_p.")



##############tabla AICc###################
result <- aictab(cand.set = Cand.models, modnames = mod.names,
                 second.ord = TRUE)

result

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt     LL
#psi.trans500_elev_p.  4 129.98       0.00   0.47   0.47 -60.63
#psi.nat500_elev_p.    4 132.38       2.40   0.14   0.61 -61.82
#psi.trans500_temp_p.  4 133.15       3.17   0.10   0.71 -62.21
#psi.trans500_lcond_p. 4 133.16       3.18   0.10   0.80 -62.22
#psi.elev_p.           3 135.01       5.03   0.04   0.84 -64.29
#psi.lcond_OD_p.       4 135.79       5.81   0.03   0.87 -63.53
#psi.nat500_temp_p.    4 135.81       5.83   0.03   0.89 -63.54
#psi.nat500_lcond_p.   4 136.28       6.30   0.02   0.91 -63.78
#psi.pH_lcond_p.       4 136.40       6.42   0.02   0.93 -63.84
#psi.elev_OD_p.        4 136.58       6.60   0.02   0.95 -63.93
#psi.trans500_pH_p.    4 137.16       7.18   0.01   0.96 -64.22
#psi.pH_p.             3 137.83       7.85   0.01   0.97 -65.70
#psi.temp_p.           3 137.92       7.94   0.01   0.98 -65.75
#psi.lcond_p.          3 138.09       8.11   0.01   0.99 -65.83
#psi.temp_OD_p.        4 138.33       8.35   0.01   0.99 -64.80
#psi.nat500_pH_p.      4 139.03       9.05   0.01   1.00 -65.15
#psi.trans500_OD_p.    4 143.98      14.00   0.00   1.00 -67.63
#psi.OD_p.             3 144.34      14.36   0.00   1.00 -68.96
#psi.nat500_OD_p.      4 145.94      15.96   0.00   1.00 -68.61
#psi.elev_2_p.         3 146.46      16.48   0.00   1.00 -70.02
#psi.trans500_p.       3 146.48      16.50   0.00   1.00 -70.02
#nulo                  2 148.49      18.51   0.00   1.00 -72.14
#psi.nat500_p.         3 149.09      19.11   0.00   1.00 -71.33

par(mfrow = c(1,1))

#####verificamos la bonda de ajuste para modelos con p constante####

obs.boot <- mb.gof.test(psi.trans500_elev_p., nsim = 1000)

print(obs.boot, digits.vals = 2, digits.chisq = 2)

#MacKenzie and Bailey goodness-of-fit for single-season occupancy model

#Pearson chi-square table:
  
#  Cohort Observed Expected Chi-square
#000       0       37    36.08       0.02
#001       0        2     1.34       0.32
#011       0        2     3.02       0.35
#100       0        2     1.34       0.32
#101       0        2     3.02       0.35
#110       0        3     3.02       0.00
#111       0        8     6.82       0.21
#NA00      1        3     3.52       0.08
#NA10      1        1     0.11       7.00

#Chi-square statistic = 10.35 
#Number of bootstrap samples = 1000
#P-value = 0.248
#
#Quantiles of bootstrapped statistics:
#  0%   25%   50%   75%  100% 
#0.29  4.13  6.65 10.27 89.42 

#Estimate of c-hat = 1.29 


result_b <- aictab(cand.set = Cand.models, modnames = mod.names,
                   second.ord = TRUE, c.hat = 1.29)

result_b

#Model selection based on QAICc:
#(c-hat estimate = 1.29)

#K  QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
#psi.trans500_elev_p.  5 105.11        0.00    0.34   0.34   -47.00
#psi.nat500_elev_p.    5 106.96        1.86    0.13   0.47   -47.93
#psi.trans500_temp_p.  5 107.56        2.46    0.10   0.57   -48.23
#psi.trans500_lcond_p. 5 107.57        2.46    0.10   0.67   -48.23
#psi.elev_p.           4 108.41        3.30    0.06   0.73   -49.84
#psi.lcond_OD_p.       5 109.61        4.51    0.04   0.77   -49.25
#psi.nat500_temp_p.    5 109.62        4.52    0.04   0.80   -49.26
#psi.nat500_lcond_p.   5 109.99        4.88    0.03   0.83   -49.44
#psi.pH_lcond_p.       5 110.08        4.98    0.03   0.86   -49.48
#psi.elev_OD_p.        5 110.22        5.12    0.03   0.89   -49.56
#psi.pH_p.             4 110.59        5.49    0.02   0.91   -50.93
#psi.temp_p.           4 110.66        5.55    0.02   0.93   -50.97
#psi.trans500_pH_p.    5 110.67        5.57    0.02   0.95   -49.78
#psi.lcond_p.          4 110.79        5.68    0.02   0.97   -51.03
#psi.temp_OD_p.        5 111.58        6.47    0.01   0.98   -50.23
#psi.nat500_pH_p.      5 112.12        7.02    0.01   0.99   -50.51
#psi.OD_p.             4 115.64       10.53    0.00   1.00   -53.46
#psi.trans500_OD_p.    5 115.96       10.86    0.00   1.00   -52.43
#psi.elev_2_p.         4 117.28       12.17    0.00   1.00   -54.28
#psi.trans500_p.       4 117.29       12.18    0.00   1.00   -54.28
#psi.nat500_OD_p.      5 117.48       12.37    0.00   1.00   -53.18
#nulo                  3 118.27       13.17    0.00   1.00   -55.92
#psi.nat500_p.         4 119.32       14.21    0.00   1.00   -55.30

##############Promedio ponderado de variables predictoras########################

modavg(cand.set = Cand.models, parm = "(Intercept)", modnames = mod.names,
       second.ord = TRUE, parm.type = "detect", c.hat = 1.29)

#Model-averaged estimate: 0.85 
#Unconditional SE: 0.38 
#95% Unconditional confidence interval: 0.11, 1.58 ##predictor confiable



modavg(cand.set = Cand.models, modnames = mod.names, second.ord = TRUE,
       parm = "trans", parm.type = "psi", c.hat = 1.29)

#Model-averaged estimate: 0.92 
#Unconditional SE: 0.49 
#95% Unconditional confidence interval: -0.04, 1.89 ##predictor no confiable

modavg(cand.set = Cand.models, modnames = mod.names, second.ord = TRUE,
       parm = "elev", parm.type = "psi", c.hat = 1.29)

#Model-averaged estimate: 1.96 
#Unconditional SE: 0.84 
#95% Unconditional confidence interval: 0.31, 3.61 

modavg(cand.set = Cand.models, modnames = mod.names, second.ord = TRUE,
       parm = "nat_500", parm.type = "psi", c.hat = 1.29)

#Model-averaged estimate: -0.71 
#Unconditional SE: 0.43 
#95% Unconditional confidence interval: -1.54, 0.13 

modavg(cand.set = Cand.models, parm = "(Intercept)", modnames = mod.names,
       second.ord = TRUE, parm.type = "psi", c.hat = 1.29)

#Model-averaged estimate: -1.05 
#Unconditional SE: 0.48 
#95% Unconditional confidence interval: -1.99, -0.12

modavg(cand.set = Cand.models, parm = "(Intercept)", modnames = mod.names,
       second.ord = TRUE, parm.type = "detect", c.hat = 1.29)

#Model-averaged estimate: 0.85 
#Unconditional SE: 0.38 
#95% Unconditional confidence interval: 0.11, 1.58

#######graficas de los valores predichos##########################
########proporción transformado a 500m #########

modavg.trans_500 <- list(psi.trans500_lcond_p., psi.trans500_p. ,
                        psi.trans500_elev_p. , psi.trans500_pH_p., 
                        psi.trans500_temp_p., psi.trans500_OD_p.)

nameavg.trans_500 <- c("psi.trans500_lcond_p.", "psi.trans500_p." ,
                      "psi.trans500_elev_p." , "psi.trans500_pH_p.", 
                      "psi.trans500_temp_p.", "psi.trans500_OD_p.")



newtrans <- data.frame(trans = seq(min(scale(amby$trans)),
                                         max(scale(amby$trans)), 
                                         length=100),
                          mean_ph=0, elev = 0, lcond=0,
                          mean_temp=0, mean_OD = 0, trans = 0)

predtrans <- modavgPred(cand.set = modavg.trans_500, modnames = nameavg.trans_500,
                           newdata = newtrans, parm.type = "psi",
                           uncond.se = "revised",type = "response", 
                           conf.level = 0.95, second.ord = "TRUE",
                           c.hat = 1.29)


##secuencia de datos de terreno transformado
trans.graf <- seq(min(amby$trans), max(amby$trans), length=100)


##extraemos los valores de los IC

avgtrans<- (predtrans$mod.avg.pred)

matxA <- predtrans$matrix.output

alower <- matxA[,"lower.CL"]

aupper <- matxA[,"upper.CL"]

####graficamos los valores predichos de terreno transformado a 500m####


par(las=1, mar=c(5,6,2,2), mfrow = c(1,3))

plot(trans.graf, avgtrans, type="l", ylim = c(0, 0.963), 
     cex.axis = 1.3, cex.lab = 1.5,  xlim = c(0, 0.963), family = "serif",
     xlab="Proporción terreno transformado 500m", las = 1,
     ylab=expression(paste("Probabilidad de ocupación ( ", 
                           italic(hat(bar(psi))),")", sep="")), bty="l")

text(0.5, 0.9, "A)", family = "serif", cex = 2.3)

lines(trans.graf, aupper, lty=2)
lines(trans.graf, alower, lty=2)

####################valores predichos de vegetación nativa################

modavg.nat <- list(psi.nat500_p., psi.nat500_temp_p. ,
                    psi.nat500_lcond_p., psi.nat500_pH_p. ,
                    psi.nat500_elev_p., psi.nat500_OD_p.)

nameavg.nat <- c("psi.nat500_p.", "psi.nat500_temp_p." ,
                  "psi.nat500_lcond_p.", "psi.nat500_pH_p." ,
                  "psi.nat500_elev_p.", "psi.nat500_OD_p.")



newnat <- data.frame(nat_500 = seq(min(scale(amby$nat_500)),
                                 max(scale(amby$nat_500)), 
                                 length=100),
                     mean_ph=0, elev = 0, lcond=0,
                     mean_temp=0, mean_OD = 0, trans = 0)

prednat <- modavgPred(cand.set = modavg.nat, modnames = nameavg.nat,
                       newdata = newnat, parm.type = "psi",
                       uncond.se = "revised",type = "response", 
                       conf.level = 0.95, second.ord = "TRUE",
                       c.hat = 1.29)


##secuencia de datos vegetación nativa
nat.graf <- seq(min(amby$nat_500), max(amby$nat_500), length=100)


##extraemos los valores de los IC

avgnat <- (prednat$mod.avg.pred)

matxB <- prednat$matrix.output

blower <- matxB[,"lower.CL"]

bupper <- matxB[,"upper.CL"]

####graficamos los valores predichos de vegetación nativa###


plot(nat.graf, avgnat, type="l", ylim = c(0, 0.963), 
     cex.axis = 1.3, cex.lab = 1.5,  xlim = c(0.03,0.963), family = "serif",
     xlab="Proporción vegetación nativa 500m", las = 1, bty = "l",
     ylab="")

text(0.5, 0.9, "B)", family = "serif", cex = 2.3)

lines(nat.graf, bupper, lty=2)
lines(nat.graf, blower, lty=2)



####################valores predichos de elevación################

modavg.elev <- list(psi.nat500_elev_p.  , psi.trans500_elev_p. ,
                    psi.elev_p.,
                    psi.elev_OD_p.)

nameavg.elev <- c("psi.nat500_elev_p."  , "psi.trans500_elev_p." ,
                   "psi.elev_p.",
                  "psi.elev_OD_p.")



newelev <- data.frame(elev = seq(min(scale(amby$elev)),
                                 max(scale(amby$elev)), 
                                 length=100),
                     nat_500 = 0, trans = 0,
                      mean_OD = 0)

predelev <- modavgPred(cand.set = modavg.elev, modnames = nameavg.elev,
                       newdata = newelev, parm.type = "psi",
                       uncond.se = "revised",type = "response", 
                       conf.level = 0.95, second.ord = "TRUE",
                       c.hat = 1.29)


##secuencia de datos de cultivo
elev.graf <- seq(min(amby$elev), max(amby$elev), length=100)


##extraemos los valores de los IC

avgelev <- (predelev$mod.avg.pred)

matxC <- predelev$matrix.output

clower <- matxC[,"lower.CL"]

cupper <- matxC[,"upper.CL"]



####graficamos los valores predichos de elevación####


plot(elev.graf, avgelev, type="l", ylim = c(0, 0.963), 
     cex.axis = 1.3, cex.lab = 1.5,  xlim = c(500,2910), family = "serif",
     xlab="Elevación (m snm)", las = 1, bty = "l",
     ylab="")

text(1800, 0.9, "C)", family = "serif", cex = 2.3)

lines(elev.graf, cupper, lty=2)
lines(elev.graf, clower, lty=2)





####Estimaciones de proporción de área ocupada####



re <- ranef(psi.trans500_elev_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)##numero de sitios

#     Estimate      2.5%     97.5%
#PAO 0.3333333 0.3333333 0.45

re <- ranef(psi.nat500_elev_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

# Estimate      2.5% 97.5%
#PAO 0.3333333 0.3333333   0.45

re <- ranef(psi.trans500_temp_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate           2.5%     97.5%
#  PAO  0.3333333 0.3333333 0.43

re <- ranef(psi.trans500_lcond_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#  Estimate      2.5% 97.5%
#PAO 0.3333333 0.3333333  0.416


re <- ranef(psi.elev_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5%     97.5%
#  PAO 0.3333333 0.3333333 0.4666667

re <- ranef(psi.lcond_OD_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5%     97.5%
#  PAO 0.3333333 0.3333333 0.46

re <- ranef(psi.nat500_temp_p. ) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5% 97.5%
#  PAO 0.3333333 0.3333333  0.45

re <- ranef(psi.nat500_lcond_p. ) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5%     97.5%
#  PAO 0.3333333 0.3333333 0.43

re <- ranef(psi.pH_lcond_p.) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5% 97.5%
#  PAO 0.3333333 0.3333333  0.45

re <- ranef(psi.elev_OD_p. ) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5%     97.5%
#  PAO 0.3333333 0.3333333 0.483

re <- ranef(psi.pH_p. ) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5%     97.5%
#  PAO 0.3333333 0.3333333 0.4666667

re <- ranef(psi.temp_p. ) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5%     97.5%
#  PAO 0.3333333 0.3333333 0.5166667

re <- ranef(psi.trans500_pH_p. ) #Random effects
EBUP <- bup(re, stat="mode")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI))/60)

#Estimate      2.5% 97.5%
#  PAO 0.3333333 0.3333333   0.4

#####promedio ponderado de la PAO####

Paoavg <- ((0.33*0.34)+(0.33*0.13)+ (0.33*0.10)+(0.33*0.10)+(0.33*0.06)+(0.33*0.04)+
             (0.33*0.04)+(0.33*0.03)+(0.33*0.03)+ (0.33*0.03)+ (0.33*0.02) +
             (0.33*0.02))
#[1] 0.3102

Paoavgup <- ((0.45*0.34)+(0.45*0.13)+ (0.43*0.10)+(0.41*0.10)+(0.46*0.06)+(0.46*0.04)+
               (0.45*0.04)+(0.43*0.03)+(0.45*0.03)+ (0.48*0.03)+ (0.51*0.02) +
               (0.4*0.02))
#[1] 0.4185


##################grafica de predicciones a lo largo de dos gradientes######


ts <- scale(seq(min(amby$trans), 
                max(amby$trans), length.out=50)) # Standardized for prediction

ep <- scale(seq(min(amby$elev), 
                max(amby$elev), length.out = 50)) # Standardised for prediction

nt <- scale(seq(min(amby$nat_500), 
                max(amby$nat_500), length.out = 50)) # Standardised for prediction



#################covariables elevación y transformado##################

pred.matrix1 <- array(NA, dim = c(50, 50))#especificación de la matriz y sus dimenciones
for(i in 1:50){
  for(j in 1:50){
    newData <- data.frame(x=0, y=0, trans=ts[i], elev=ep[j], iLength=0)##ingreso de datos especificados para las predicciones conjuntas
    pred.matrix1[i,j] <- predict(psi.trans500_elev_p.,
                                 type="state", newdata=newData)[1,1]##generación de las predicciones a partir del modelo
  }
}


par(mfrow = c(1, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)

mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

image(y=seq(0, 1, length.out = 50),
      x=seq(500, 3000, length.out = 50), family = "serif",
      z=pred.matrix1, col = mapPalette(50), axes = F,
      xlab =  "Elevación (m snm)", 
      ylab =  "Proporción de terreno transformado")

contour(y=seq(0, 1, length.out = 50),
        x=seq(500, 3000, length.out = 50),
        z=pred.matrix1,family = "serif", 
        add = T, col = "blue", labcex = 1.5,
        lwd = 2)

axis(2, at = seq(0, 1, by = 0.10), family = "serif")
axis(1, at = seq(500, 3000, by = 500), family = "serif", las = 1)
box()
points(amby$elev, amby$trans, pch=1, cex=1.5)

points(amby$elev[amby$s3=="1"], amby$nat_500[amby$s3=="1"], 
       pch=19, cex=1.5, col = "black")

points(amby$elev[amby$s2=="1"], amby$nat_500[amby$s2=="1"], 
       pch=19, cex=1.5, col = "black")

points(amby$elev[amby$s1=="1"], amby$nat_500[amby$s1=="1"], 
       pch=19, cex=1.5, col = "black")

#################covariables elevación y nativa##################

pred.matrix2 <- array(NA, dim = c(50, 50))#especificación de la matriz y sus dimenciones
for(i in 1:50){
  for(j in 1:50){
    newData <- data.frame(x=0, y=0, elev=ep[i], nat_500=nt[j], iLength=0)##ingreso de datos especificados para las predicciones conjuntas
    pred.matrix2[i,j] <- predict(psi.nat500_elev_p.,
                                 type="state", newdata=newData)[1,1]##generación de las predicciones a partir del modelo
  }
}



image(y=seq(0, 1, length.out = 50),
      x=seq(500, 3000, length.out = 50), family = "serif",
      z=pred.matrix2, col = mapPalette(50), axes = F,
      xlab =  "Elevación (m snm)", 
      ylab =  "Proporción de vegetación nativa")

contour(y=seq(0, 1, length.out = 50),
        x=seq(500, 3000, length.out = 50),
        z=pred.matrix2,family = "serif", 
        add = T, col = "blue", labcex = 1.5,
        lwd = 2)

axis(2, at = seq(0, 1, by = 0.10), family = "serif")
axis(1, at = seq(500, 3000, by = 500), family = "serif", las = 1)
box()
points(amby$elev, amby$nat_500, pch=1, cex=1.5)

points(amby$elev[amby$s3=="1"], amby$trans[amby$s3=="1"], 
       pch=19, cex=1.5)

points(amby$elev[amby$s2=="1"], amby$trans[amby$s2=="1"], 
       pch=19, cex=1.5)

points(amby$elev[amby$s1=="1"], amby$trans[amby$s1=="1"], 
       pch=19, cex=1.5)



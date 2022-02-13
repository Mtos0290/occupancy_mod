#####################################################################################
#### Analisis de abundancia de A. ordinarium 60 sitios  en febrero a mayo de 2018####
#### analisis rems fecha: julio 2019                                           ####
#####################################################################################

R
library(unmarked)
library(plotrix)
library(AICcmodavg)
library(usdm)
library(corrplot)
library(AHMbook)

load("abundancias.secas040719.RData")

###### Importacion y manipulacion de la base de datos###############

amby <- read.table("ambystoma_conteos_b.csv", sep= ",", header = T)

summary(amby)

names(amby)

############estimacion de los promedios de variables para secas#####

amby$mean_temp <- rowMeans(amby[,8:10])

amby$mean_ph <-  rowMeans(amby[,14:16])

amby$mean_cond <-  rowMeans(amby[,20:22])

amby$mean_OD <-  rowMeans(amby[,26:28])


names(amby)

###historiales de detección en secas####

amby.y <- amby[,2:4]


#####transfomacion a logaritmo de coductividad#####
hist(amby$mean_cond)

amby$lcond <- log(amby$mean_cond)

hist(amby$lcond)

names(amby)


####definición de las variables predictoras para psi en secas#####
amby.cov <- data.frame(scale(amby[c(38,45:47, 49,50,52,53)]))


summary(amby.cov)


####verificamos el grado de correlación entre variables###########

vif(amby.cov)

#   Variables      VIF
#1       elev 3.521673 el VIF es menor a 10 en todos los casos
#2 livest_500 2.516010
#3   cult_500 4.315668
#4    nat_500 7.287272
#5  mean_temp 3.372889
#6    mean_ph 2.909613
#7    mean_OD 2.163078
#8      lcond 2.529572


corrplot(cor(amby.cov), method="number", type = "upper",
         family = "serif", tl.cex = 1.5, number.cex = 1.2,
         tl.offset = T, tl.col = "black")
##elevación está correlacionada con temperatura y pH, así como 
#la vegetación nativa con los otros usos de suelo. Y la temperatura con
#conductividad


####estandarización de variables predictoras de p############
temp <- data.frame(scale(amby[,8:10]))
prof <- data.frame(scale(amby[,32:34]))
cobr <- data.frame(scale(amby[,39:41]))




amby.obs <- list(temp=temp, prof=prof, cobr=cobr)

ambyUMF <- unmarkedFramePCount(y = amby.y, siteCovs = amby.cov, obsCovs = amby.obs)

summary(ambyUMF)

plot(ambyUMF, col.regions = colorRampPalette(c('gray', 'darkred')))

#unmarkedFrame Object

#60 sites
#Maximum number of observations per site: 3 
#Mean number of observations per site: 2.93 
#Sites with at least one detection: 17 

#Tabulation of y observations:
#  0    1    2    3    4    5    6    8    9 <NA> 
#  140   14    9    3    1    3    3    2    1    4 

################variables para definir la detectabilidad ####

n._p.temp <- pcount(~temp ~1, data = ambyUMF)

n._p.prof <- pcount(~ prof ~1, data = ambyUMF)

n._p.cobr <- pcount(~cobr ~1, data = ambyUMF)

n._p. <- pcount(~1 ~1, data = ambyUMF)

(btd <- backTransform(n._p., type="det"))

confint(btd, level = 0.95)

n._p.tempprof <- pcount(~temp  + prof ~1, data = ambyUMF)

n._p.profcobr <- pcount(~ prof + cobr~1, data = ambyUMF)

n._p.cobrtemp <- pcount(~cobr + temp~1, data = ambyUMF)

n._p.tempXprof <- pcount(~temp * prof -temp - prof ~1, data = ambyUMF)

n._p.profXcobr <- pcount(~ prof * cobr -prof -cobr ~1, data = ambyUMF)

n._p.cobrXtemp <- pcount(~cobr * temp -cobr - temp ~1, data = ambyUMF)


#####selección de variable predictora para abundancia#########

Cand.models_p <- list(n._p.,  
                      n._p.temp, n._p.prof, 
                      n._p.cobr,
                      n._p.tempprof,
                      n._p.profcobr,
                      n._p.cobrtemp,
                      n._p.tempXprof,
                      n._p.profXcobr,
                      n._p.cobrXtemp)


mod.names_p <- c("n._p.", 
                 "n._p.temp", "n._p.prof", 
                 "n._p.cobr",
                 "n._p.tempprof",
                 "n._p.profcobr",
                 "n._p.cobrtemp",
                 "n._p.tempXprof",
                 "n._p.profXcobr",
                 "n._p.cobrXtemp")

aictab(cand.set = Cand.models_p, modnames = mod.names_p,
       second.ord = TRUE)

#Model selection based on AICc:
  
#                   K   AICc Delta_AICc AICcWt Cum.Wt      LL
#n._p.cobrXtemp      3 380.36       0.00   0.95   0.95 -186.97#modelos con más parsimonioso
#n._p.tempXprof      3 387.07       6.70   0.03   0.98 -190.32
#n._p.prof           3 390.51      10.15   0.01   0.99 -192.04
#n._p.profcobr       4 392.33      11.97   0.00   0.99 -191.80
#n._p.profXcobr      3 392.55      12.19   0.00   0.99 -193.06
#n._p.tempprof       4 392.69      12.33   0.00   1.00 -191.98
#n._p.               2 394.11      13.75   0.00   1.00 -194.95
#n._p.cobr           3 396.13      15.76   0.00   1.00 -194.85
#n._p.temp           3 396.22      15.85   0.00   1.00 -194.89
#n._p.cobrtemp       4 397.39      17.02   0.00   1.00 -194.33


###########modelos N-mix##################


########################vegetación nativa######################

n.nat_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~nat_500, data = ambyUMF)

n.nat_temp_p. <- pcount(~1 ~nat_500 + mean_temp, data = ambyUMF)

n.nat_lcond_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~nat_500 + lcond, data = ambyUMF)

n.nat_pH_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~nat_500 + mean_ph, data = ambyUMF)

n.nat_OD_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~nat_500 + mean_OD, data = ambyUMF)

n.nat_elev_p.cobXtemp <- pcount(~1 ~nat_500 + elev, data = ambyUMF)

########################potrero##########################
n.lives_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~livest_500, data = ambyUMF)

n.lives_temp_p. <- pcount(~1 ~livest_500 + mean_temp, data = ambyUMF)

n.lives_lcond_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~livest_500 + lcond, data = ambyUMF)

n.lives_pH_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~livest_500 + mean_ph, data = ambyUMF)

n.lives_OD_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~livest_500 + mean_OD, data = ambyUMF)

n.lives_elev_p.cobXtemp <- pcount(~1 ~livest_500 + elev, data = ambyUMF)

###############################Cultivos###########################

n.cult_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~cult_500, data = ambyUMF)

n.cult_temp_p. <- pcount(~1 ~cult_500 + mean_temp, data = ambyUMF)

n.cult_lcond_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~cult_500 + lcond, data = ambyUMF)

n.cult_pH_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~cult_500 + mean_ph, data = ambyUMF)

n.cult_OD_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~cult_500 + mean_OD, data = ambyUMF)

n.cult_elev_p.cobXtemp <- pcount(~1 ~cult_500 + elev, data = ambyUMF)

#########################variables locales########################

n.temp_p.cobXtemp <- pcount(~1 ~ mean_temp, data=ambyUMF)

n.lcond_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~ lcond, data = ambyUMF)

n.pH_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~ mean_ph, data = ambyUMF)

n.OD_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~ mean_OD, data = ambyUMF)

n.elev_p.cobXtemp <- pcount(~1 ~ elev, data = ambyUMF)

n.temp_OD_p.cobXtemp <- pcount(~1 ~ mean_temp + mean_OD, data = ambyUMF)

n.lcond_OD_p.cobXtemp <- pcount(~cobr * temp -cobr -temp ~ lcond + mean_OD, data = ambyUMF)

NULO <- pcount(~1 ~1, data = ambyUMF)


n.elev_2_p.cobXtemp <- pcount(~1 ~ I(elev^2), data = ambyUMF)

n.lcond_elev_p.cobXtemp <- pcount(~1 ~ lcond + elev, data = ambyUMF)

n.elev_OD_p.cobXtemp <- pcount(~1 ~ elev + mean_OD, data = ambyUMF)

n.temp_pH_p.cobXtemp <- pcount(~1 ~ mean_temp + mean_ph, data = ambyUMF)

n.pH_OD_p.cobXtemp <- pcount(~cobr:temp ~ mean_ph + mean_OD, data = ambyUMF)

n.pH_lcond_p.cobXtemp <- pcount(~cobr:temp ~ mean_ph + lcond, data = ambyUMF)

#################selección modelos #######
##############lista de modelos############


Cand.models_n <- list(n.nat_p.cobXtemp, n.nat_temp_p., n.nat_lcond_p.cobXtemp,
                      n.nat_pH_p.cobXtemp,
                      n.nat_OD_p.cobXtemp,n.nat_elev_p.cobXtemp,
                      n.lives_p.cobXtemp, n.lives_temp_p., 
                      n.lives_lcond_p.cobXtemp, 
                      n.lives_pH_p.cobXtemp, n.lives_OD_p.cobXtemp, 
                      n.lives_elev_p.cobXtemp, n.cult_p.cobXtemp, 
                      n.cult_temp_p., n.cult_lcond_p.cobXtemp,
                      n.cult_pH_p.cobXtemp,n.cult_OD_p.cobXtemp,
                      n.cult_elev_p.cobXtemp,
                      n.temp_p.cobXtemp, n.lcond_p.cobXtemp,
                      n.pH_p.cobXtemp, n.OD_p.cobXtemp,
                      n.elev_p.cobXtemp, n.temp_OD_p.cobXtemp,
                      n.lcond_OD_p.cobXtemp,
                      NULO, n.elev_2_p.cobXtemp,
                      n.lcond_elev_p.cobXtemp, n.elev_OD_p.cobXtemp, 
                      n.temp_pH_p.cobXtemp, 
                      n.pH_OD_p.cobXtemp, n.pH_lcond_p.cobXtemp)

##############nombres de candidatos###########

mod.names_n <- c("n.nat_p.cobXtemp", "n.nat_temp_p.", "n.nat_lcond_p.cobXtemp",
                 "n.nat_pH_p.cobXtemp",
                 "n.nat_OD_p.cobXtemp","n.nat_elev_p.cobXtemp",
                 "n.lives_p.cobXtemp", "n.lives_temp_p.", 
                 "n.lives_lcond_p.cobXtemp", 
                 "n.lives_pH_p.cobXtemp", "n.lives_OD_p.cobXtemp", 
                 "n.lives_elev_p.cobXtemp", "n.cult_p.cobXtemp", 
                 "n.cult_temp_p.", "n.cult_lcond_p.cobXtemp",
                 "n.cult_pH_p.cobXtemp","n.cult_OD_p.cobXtemp",
                 "n.cult_elev_p.cobXtemp",
                 "n.temp_p.cobXtemp", "n.lcond_p.cobXtemp",
                 "n.pH_p.cobXtemp", "n.OD_p.cobXtemp",
                 "n.elev_p.cobXtemp", "n.temp_OD_p.cobXtemp",
                 "n.lcond_OD_p.cobXtemp",
                 "NULO", "n.elev_2_p.cobXtemp",
                 "n.lcond_elev_p.cobXtemp", "n.elev_OD_p.cobXtemp", 
                 "n.temp_pH_p.cobXtemp", 
                 "n.pH_OD_p.cobXtemp", "n.pH_lcond_p.cobXtemp")



##############tabla AICc###################
result.n <- aictab(cand.set = Cand.models_n, modnames = mod.names_n,
                 second.ord = TRUE)

result.n

#Model selection based on AICc:

#K   AICc Delta_AICc AICcWt Cum.Wt      LL
#n.cult_elev_p.cobXtemp   4 310.59       0.00   0.99   0.99 -150.93
#n.cult_pH_p.cobXtemp     5 320.85      10.26   0.01   1.00 -154.87
#n.lcond_OD_p.cobXtemp    5 325.21      14.62   0.00   1.00 -157.05
#n.cult_lcond_p.cobXtemp  5 326.08      15.50   0.00   1.00 -157.49
#n.pH_lcond_p.cobXtemp    5 328.75      18.16   0.00   1.00 -158.82
#n.lcond_elev_p.cobXtemp  4 329.69      19.10   0.00   1.00 -160.48
#n.pH_OD_p.cobXtemp       5 329.80      19.22   0.00   1.00 -159.35
#n.nat_elev_p.cobXtemp    4 330.40      19.81   0.00   1.00 -160.83
#n.pH_p.cobXtemp          4 332.57      21.98   0.00   1.00 -161.92
#n.nat_pH_p.cobXtemp      5 332.69      22.10   0.00   1.00 -160.79
#n.lives_pH_p.cobXtemp    5 334.32      23.73   0.00   1.00 -161.60
#n.elev_OD_p.cobXtemp     4 335.81      25.22   0.00   1.00 -163.54
#n.lives_elev_p.cobXtemp  4 337.89      27.30   0.00   1.00 -164.58
#n.elev_p.cobXtemp        3 339.53      28.95   0.00   1.00 -166.55
#n.cult_OD_p.cobXtemp     5 340.03      29.44   0.00   1.00 -164.46
#n.nat_lcond_p.cobXtemp   5 340.04      29.45   0.00   1.00 -164.46
#n.lives_lcond_p.cobXtemp 5 343.28      32.70   0.00   1.00 -166.09
#n.lcond_p.cobXtemp       4 344.09      33.50   0.00   1.00 -167.68
#n.temp_pH_p.cobXtemp     4 345.19      34.61   0.00   1.00 -168.23
#n.OD_p.cobXtemp          4 350.49      39.90   0.00   1.00 -170.88
#n.nat_OD_p.cobXtemp      5 352.05      41.46   0.00   1.00 -170.47
#n.lives_OD_p.cobXtemp    5 352.10      41.51   0.00   1.00 -170.49
#n.cult_temp_p.           4 356.06      45.47   0.00   1.00 -173.66
#n.cult_p.cobXtemp        4 357.83      47.24   0.00   1.00 -174.55
#n.temp_OD_p.cobXtemp     4 362.95      52.36   0.00   1.00 -177.11
#n.nat_temp_p.            4 371.27      60.68   0.00   1.00 -181.27
#n.lives_temp_p.          4 374.07      63.48   0.00   1.00 -182.67
#n.temp_p.cobXtemp        3 375.86      65.27   0.00   1.00 -184.72
#n.nat_p.cobXtemp         4 379.26      68.67   0.00   1.00 -185.27
#n.lives_p.cobXtemp       4 381.18      70.60   0.00   1.00 -186.23
#NULO                     2 394.11      83.53   0.00   1.00 -194.95
#n.elev_2_p.cobXtemp      3 394.33      83.74   0.00   1.00 -193.95


summary(n.cult_elev_p.cobXtemp)

#Call:
#pcount(formula = ~1 ~ cult_500 + elev, data = ambyUMF)

#Abundance (log-scale):
#  Estimate    SE     z  P(>|z|)
#(Intercept)   -0.685 0.297 -2.31 2.10e-02
#cult_500       0.624 0.110  5.67 1.43e-08
#elev           1.666 0.256  6.51 7.60e-11

#Detection (logit-scale):
#  Estimate    SE     z P(>|z|)
#-0.44 0.279 -1.58   0.115

#AIC: 309.8602 
#Number of sites: 60
#optim convergence code: 0
#optim iterations: 46 
#Bootstrap iterations: 0

####verificamos de la bondad de ajuste del modelo mas parsimonioso####

obs.bootP <- Nmix.gof.test(n.cult_elev_p.cobXtemp, nsim = 1000)##20 minutos

print(obs.bootP, digits.vals = 2, digits.chisq = 2)
#Estimate of c-hat = 2.37

result.Q <- aictab(cand.set = Cand.models_n, modnames = mod.names_n,
                   second.ord = TRUE, c.hat=2.37)

result.Q

##Model selection based on QAICc:
#(c-hat estimate = 2.37)

#K  QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
#n.cult_elev_p.cobXtemp   5 138.48        0.00    0.85   0.85   -63.68
#n.cult_pH_p.cobXtemp     6 144.28        5.80    0.05   0.89   -65.35
#n.lcond_OD_p.cobXtemp    6 146.11        7.64    0.02   0.91   -66.26
#n.cult_lcond_p.cobXtemp  6 146.48        8.01    0.02   0.93   -66.45
#n.lcond_elev_p.cobXtemp  5 146.54        8.06    0.02   0.94   -67.71
#n.nat_elev_p.cobXtemp    5 146.84        8.36    0.01   0.96   -67.86
#n.pH_lcond_p.cobXtemp    6 147.61        9.13    0.01   0.96   -67.01
#n.pH_p.cobXtemp          5 147.75        9.28    0.01   0.97   -68.32
#n.pH_OD_p.cobXtemp       6 148.05        9.58    0.01   0.98   -67.23
#n.elev_OD_p.cobXtemp     5 149.12       10.64    0.00   0.98   -69.00
#n.nat_pH_p.cobXtemp      6 149.27       10.79    0.00   0.99   -67.84
#n.elev_p.cobXtemp        4 149.28       10.80    0.00   0.99   -70.28
#n.lives_pH_p.cobXtemp    6 149.96       11.48    0.00   0.99   -68.19
#n.lives_elev_p.cobXtemp  5 150.00       11.52    0.00   1.00   -69.44
#n.cult_OD_p.cobXtemp     6 152.37       13.89    0.00   1.00   -69.39
#n.nat_lcond_p.cobXtemp   6 152.37       13.89    0.00   1.00   -69.39
#n.lcond_p.cobXtemp       5 152.61       14.14    0.00   1.00   -70.75
#n.temp_pH_p.cobXtemp     5 153.08       14.60    0.00   1.00   -70.98
#n.lives_lcond_p.cobXtemp 6 153.74       15.26    0.00   1.00   -70.08
#n.OD_p.cobXtemp          5 155.32       16.84    0.00   1.00   -72.10
#n.nat_OD_p.cobXtemp      6 157.44       18.96    0.00   1.00   -71.93
#n.lives_OD_p.cobXtemp    6 157.46       18.98    0.00   1.00   -71.94
#n.cult_temp_p.           5 157.66       19.18    0.00   1.00   -73.28
#n.cult_p.cobXtemp        5 158.41       19.93    0.00   1.00   -73.65
#n.temp_OD_p.cobXtemp     5 160.57       22.09    0.00   1.00   -74.73
#n.nat_temp_p.            5 164.08       25.60    0.00   1.00   -76.48
#n.temp_p.cobXtemp        4 164.61       26.13    0.00   1.00   -77.94
#n.lives_temp_p.          5 165.26       26.79    0.00   1.00   -77.08
#n.nat_p.cobXtemp         5 167.45       28.98    0.00   1.00   -78.17
#n.lives_p.cobXtemp       5 168.27       29.79    0.00   1.00   -78.58
#NULO                     3 170.94       32.47    0.00   1.00   -82.26
#n.elev_2_p.cobXtemp      4 172.40       33.92    0.00   1.00   -81.84

n.cult_elev_p.cobXtempNB <- pcount(~1~cult_500 + elev,
                                   data = ambyUMF, mixture = "NB")


n.cult_elev_p.cobXtempZIP <- pcount(~1~cult_500 + elev,
                                   data = ambyUMF, mixture = "ZIP")

###tabla de AICc
Cand.models2 <- list(n.cult_elev_p.cobXtemp, n.cult_elev_p.cobXtempNB,
                     n.cult_elev_p.cobXtempZIP)

mod.names2 <- c("n.cult_elev_p.cobXtemp", "n.cult_elev_p.cobXtempNB",
                "n.cult_elev_p.cobXtempZIP")

aictab(cand.set = Cand.models2, modnames = mod.names2, second.ord = T)

#Model selection based on AICc:

#K   AICc Delta_AICc AICcWt Cum.Wt      LL
#n.cult_elev_p.cobXtempNB  6 250.46       0.00   0.77   0.77 -118.44
#n.cult_elev_p.cobXtempZIP 6 252.88       2.42   0.23   1.00 -119.65
#n.cult_elev_p.cobXtemp    5 300.27      49.81   0.00   1.00 -144.58

#Estimamos la bondad de ajuste de los modelos con distribución NB y ZIP

(gof.NB <- Nmix.gof.test(n.cult_elev_p.cobXtempNB, nsim = 1000))
#Estimate of c-hat =0.89

(gof.ZIP <- Nmix.gof.test(n.cult_elev_p.cobXtempZIP, nsim = 1000))
#Estimate of c-hat = 1.34

save.image("abundancias.secas040719.RData")

#revisamos el ajuste de los residuales
plot_Nmix_resi(n.cult_elev_p.cobXtemp, n.cult_elev_p.cobXtempNB, n.cult_elev_p.cobXtempZIP)

#el  modelo con distribución NB tiene los valores ajustados en menor distribución que 
#los observados en comparación con los modelos poisson y ZIP

###estimación de la desviación de la raíz cuadrada media###

(RP <- sqrt(mean(as.matrix(((amby.y) - fitted(n.cult_elev_p.cobXtemp))^2), na.rm = T)))#Poisson
(RNB <- sqrt(mean(as.matrix(((amby.y) - fitted(n.cult_elev_p.cobXtempNB))^2), na.rm = T)))#NB
(RZIP <- sqrt(mean(as.matrix(((amby.y) - fitted(n.cult_elev_p.cobXtempZIP))^2), na.rm = T)))#ZIP
##hay que agregar as.matrix 
#[1] 1.389648
#[1] 1.418416
#[1] 1.513524
#


summary(n.cult_elev_p.cobXtempZIP)
#Call:
#  pcount(formula = ~1 ~ cult_500 + elev, data = ambyUMF, mixture = "ZIP")

#Abundance (log-scale):
#  Estimate    SE    z  P(>|z|)
#(Intercept)    0.534 0.712 0.75 4.53e-01
#cult_500       0.428 0.105 4.07 4.60e-05
#elev           2.229 0.340 6.56 5.45e-11

#Detection (logit-scale):
#  Estimate    SE     z P(>|z|)
#-1.65 0.668 -2.46  0.0137

#Zero-inflation (logit-scale):
#  Estimate    SE       z P(>|z|)
#-0.0275 0.409 -0.0672   0.946

#AIC: 252.6887 
#Number of sites: 60
#optim convergence code: 0
#optim iterations: 25 
#Bootstrap iterations: 0 

#probamos modelos con sólo una variable para la detección



n.cult_elev_p.tempZIP <- pcount(~1 ~cult_500 + elev, data = ambyUMF,
                                mixture = "ZIP")#modelo solo con temperatura para p

summary(n.cult_elev_p.cobZIP)



###############selección de modelos con distribución ZIP #####
########################vegetación nativa######################

n.nat_p.tempZIP <- pcount(~temp ~nat_500, data = ambyUMF, mixture = "ZIP")

n.nat_temp_p.tempZIP <- pcount(~1 ~nat_500 + mean_temp, data = ambyUMF,
                         mixture = "ZIP")

n.nat_lcond_p.tempZIP <- pcount(~temp ~nat_500 + lcond, data = ambyUMF,
                          mixture = "ZIP")

n.nat_pH_p.tempZIP <- pcount(~temp ~nat_500 + mean_ph, data = ambyUMF,
                       mixture = "ZIP")

n.nat_OD_p.tempZIP <- pcount(~temp ~nat_500 + mean_OD, data = ambyUMF,
                       mixture = "ZIP")

n.nat_elev_p.tempZIP <- pcount(~1 ~nat_500 + elev, data = ambyUMF,
                         mixture = "ZIP")

########################potrero##########################
n.lives_p.tempZIP <- pcount(~temp ~livest_500, data = ambyUMF,
                      mixture = "ZIP")

n.lives_temp_p.tempZIP <- pcount(~1 ~livest_500 + mean_temp, data = ambyUMF,
                           mixture = "ZIP")

n.lives_lcond_p.tempZIP <- pcount(~temp ~livest_500 + lcond, data = ambyUMF,
                            mixture = "ZIP")

n.lives_pH_p.tempZIP <- pcount(~temp ~livest_500 + mean_ph, data = ambyUMF,
                         mixture = "ZIP")

n.lives_OD_p.tempZIP <- pcount(~temp ~livest_500 + mean_OD, data = ambyUMF,
                         mixture = "ZIP")

n.lives_elev_p.tempZIP <- pcount(~1 ~livest_500 + elev, data = ambyUMF,
                           mixture = "ZIP")

###############################Cultivos###########################

n.cult_p.tempZIP <- pcount(~temp ~cult_500, data = ambyUMF,
                     mixture = "ZIP")

n.cult_temp_p.tempZIP <- pcount(~1 ~cult_500 + mean_temp, data = ambyUMF,
                          mixture = "ZIP")

n.cult_lcond_p.tempZIP <- pcount(~temp ~cult_500 + lcond, data = ambyUMF,
                           mixture = "ZIP")

n.cult_pH_p.tempZIP <- pcount(~temp ~cult_500 + mean_ph, data = ambyUMF,
                        mixture = "ZIP")

n.cult_OD_p.tempZIP <- pcount(~temp ~cult_500 + mean_OD, data = ambyUMF,
                        mixture = "ZIP")



#########################variables locales########################

n.temp_p.tempZIP <- pcount(~1 ~ mean_temp, data=ambyUMF,
                     mixture = "ZIP")

n.lcond_p.tempZIP <- pcount(~temp ~ lcond, data = ambyUMF,
                      mixture = "ZIP")

n.pH_p.tempZIP <- pcount(~temp ~ mean_ph, data = ambyUMF,
                   mixture = "ZIP")

n.OD_p.tempZIP <- pcount(~temp ~ mean_OD, data = ambyUMF,
                   mixture = "ZIP")

n.elev_p.tempZIP <- pcount(~1 ~ elev, data = ambyUMF,
                     mixture = "ZIP")

n.temp_OD_p.tempZIP <- pcount(~1 ~ mean_temp + mean_OD, data = ambyUMF,
                        mixture = "ZIP")

n.lcond_OD_p.tempZIP <- pcount(~temp ~ lcond + mean_OD, data = ambyUMF,
                         mixture = "ZIP")

NULO.ZIP <- pcount(~1 ~1, data = ambyUMF,
               mixture = "ZIP")

n.cult_elev_p.ZIP <- pcount(~1 ~elev + cult_500, data = ambyUMF,
          mixture = "ZIP")


n.elev_2_p.tempZIP <- pcount(~temp ~ I(elev^2), data = ambyUMF,
                       mixture = "ZIP")

n.lcond_elev_p.tempZIP <- pcount(~1 ~ lcond + elev, data = ambyUMF,
                           mixture = "ZIP")

n.elev_OD_p.tempZIP <- pcount(~1 ~ elev + lcond, data = ambyUMF,
                        mixture = "ZIP")

n.temp_pH_p.tempZIP <- pcount(~1 ~ mean_temp + mean_ph, data = ambyUMF,
                        mixture = "ZIP")

n.pH_OD_p.tempZIP <- pcount(~temp ~ mean_ph + mean_OD, data = ambyUMF,
                      mixture = "ZIP")

n.pH_lcond_p.tempZIP <- pcount(~temp ~ mean_ph + lcond, data = ambyUMF,
                         mixture = "ZIP")

#################selección modelos #######
##############lista de modelos############


Cand.models.ZIP <- list(n.nat_p.tempZIP, n.nat_temp_p.tempZIP, 
                      n.nat_lcond_p.tempZIP, n.nat_pH_p.tempZIP,
                      n.nat_OD_p.tempZIP, n.nat_elev_p.tempZIP, 
                      n.lives_p.tempZIP, n.lives_temp_p.tempZIP, 
                      n.lives_lcond_p.tempZIP, n.lives_pH_p.tempZIP,  
                      n.lives_OD_p.tempZIP,  n.lives_elev_p.tempZIP,  
                      n.cult_p.tempZIP, n.cult_temp_p.tempZIP, 
                      n.cult_lcond_p.tempZIP, n.cult_pH_p.tempZIP,
                      n.cult_OD_p.tempZIP, 
                      n.temp_p.tempZIP, 
                      n.lcond_p.tempZIP, 
                      n.pH_p.tempZIP, 
                      n.OD_p.tempZIP, 
                      n.elev_p.tempZIP, 
                      n.temp_OD_p.tempZIP,
                      n.lcond_OD_p.tempZIP,
                      NULO.ZIP, n.elev_2_p.tempZIP,n.lcond_elev_p.tempZIP,
                      n.elev_OD_p.tempZIP, n.temp_pH_p.tempZIP, 
                      n.pH_OD_p.tempZIP, n.pH_lcond_p.tempZIP, n.cult_elev_p.ZIP)

##############nombres de candidatos###########

mod.names.ZIP <- c("n.nat_p.tempZIP", "n.nat_temp_p.tempZIP", 
                "n.nat_lcond_p.tempZIP", "n.nat_pH_p.tempZIP",
                "n.nat_OD_p.tempZIP", "n.nat_elev_p.tempZIP", 
                "n.lives_p.tempZIP", "n.lives_temp_p.tempZIP", 
                "n.lives_lcond_p.tempZIP", "n.lives_pH_p.tempZIP",  
                "n.lives_OD_p.tempZIP",  "n.lives_elev_p.tempZIP",  
                "n.cult_p.tempZIP", "n.cult_temp_p.tempZIP", 
                "n.cult_lcond_p.tempZIP", "n.cult_pH_p.tempZIP",
                "n.cult_OD_p.tempZIP", 
                "n.temp_p.tempZIP", 
                "n.lcond_p.tempZIP", 
                "n.pH_p.tempZIP", 
                "n.OD_p.tempZIP", 
                "n.elev_p.tempZIP", 
                "n.temp_OD_p.tempZIP",
                "n.lcond_OD_p.tempZIP",
                "NULO", "n.elev_2_p.tempZIP",
                "n.lcond_elev_p.tempZIP",
                "n.elev_OD_p.tempZIP", "n.temp_pH_p.tempZIP", 
                "n.pH_OD_p.tempZIP", "n.pH_lcond_p.tempZIP", "n.cult_elev_p.ZIP")





result.ZIP <- aictab(cand.set = Cand.models.ZIP, modnames = mod.names.ZIP,
                   second.ord = TRUE)

result.ZIP



####
#

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt      LL
#n.cult_elev_p.ZIP       5 253.80       0.00   0.94   0.94 -121.34
#n.nat_elev_p.tempZIP    5 259.56       5.76   0.05   0.99 -124.23
#n.lives_elev_p.tempZIP  5 264.78      10.98   0.00   1.00 -126.84
#n.elev_p.tempZIP        4 267.74      13.94   0.00   1.00 -129.51
#n.lcond_elev_p.tempZIP  5 269.81      16.01   0.00   1.00 -129.35
#n.elev_OD_p.tempZIP     5 269.81      16.01   0.00   1.00 -129.35
#n.lcond_OD_p.tempZIP    6 273.89      20.09   0.00   1.00 -130.15
#n.elev_2_p.tempZIP      5 274.98      21.18   0.00   1.00 -131.94
#n.lcond_p.tempZIP       5 286.94      33.14   0.00   1.00 -137.92
#n.lives_lcond_p.tempZIP 6 287.17      33.37   0.00   1.00 -136.79
#n.pH_lcond_p.tempZIP    6 288.65      34.85   0.00   1.00 -137.53
#n.pH_p.tempZIP          5 289.24      35.44   0.00   1.00 -139.06
#n.nat_lcond_p.tempZIP   6 289.41      35.61   0.00   1.00 -137.91
#n.cult_lcond_p.tempZIP  6 289.42      35.62   0.00   1.00 -137.92
#n.pH_OD_p.tempZIP       6 290.02      36.22   0.00   1.00 -138.22
#n.temp_pH_p.tempZIP     5 290.91      37.11   0.00   1.00 -139.90
#n.lives_pH_p.tempZIP    6 291.39      37.59   0.00   1.00 -138.90
#n.nat_pH_p.tempZIP      6 291.53      37.73   0.00   1.00 -138.97
#n.cult_pH_p.tempZIP     6 291.58      37.78   0.00   1.00 -139.00
#n.temp_p.tempZIP        4 292.35      38.55   0.00   1.00 -141.81
#n.lives_p.tempZIP       5 292.76      38.96   0.00   1.00 -140.83
#n.OD_p.tempZIP          5 292.77      38.97   0.00   1.00 -140.83
#n.nat_p.tempZIP         5 292.98      39.18   0.00   1.00 -140.94
#n.cult_p.tempZIP        5 292.99      39.19   0.00   1.00 -140.94
#NULO                    3 293.50      39.70   0.00   1.00 -143.54
#n.nat_temp_p.tempZIP    5 293.93      40.13   0.00   1.00 -141.41
#n.cult_temp_p.tempZIP   5 294.14      40.34   0.00   1.00 -141.51
#n.temp_OD_p.tempZIP     5 294.24      40.44   0.00   1.00 -141.56
#n.lives_temp_p.tempZIP  5 294.36      40.56   0.00   1.00 -141.63
#n.nat_OD_p.tempZIP      6 294.60      40.80   0.00   1.00 -140.51
#n.cult_OD_p.tempZIP     6 294.63      40.83   0.00   1.00 -140.52
#n.lives_OD_p.tempZIP    6 294.82      41.02   0.00   1.00 -140.62##





print(gof.ZIP, digits.vals = 2, digits.chisq = 2)
#Estimate of c-hat = 1.34 hay un poco de sobredispersión

save.image("abundancias.secas040719.RData")

result.ZIPQ <- aictab(cand.set = Cand.models.ZIP, modnames = mod.names.ZIP,
                     second.ord = TRUE, c.hat =1.34)

result.ZIPQ


#Model selection based on QAICc:
#(c-hat estimate = 1.34)

#K  QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
#n.cult_elev_p.ZIP       6 194.70        0.00    0.87   0.87   -90.56
#n.nat_elev_p.tempZIP    6 199.00        4.30    0.10   0.97   -92.71
#n.lives_elev_p.tempZIP  6 202.89        8.20    0.01   0.99   -94.65
#n.elev_p.tempZIP        5 204.40        9.71    0.01   0.99   -96.65
#n.lcond_elev_p.tempZIP  6 206.64       11.95    0.00   1.00   -96.53
#n.elev_OD_p.tempZIP     6 206.64       11.95    0.00   1.00   -96.53
#n.lcond_OD_p.tempZIP    7 210.41       15.71    0.00   1.00   -97.13
#n.elev_2_p.tempZIP      6 210.50       15.81    0.00   1.00   -98.46
#n.lcond_p.tempZIP       6 219.43       24.73    0.00   1.00  -102.92
#n.lives_lcond_p.tempZIP 7 220.32       25.63    0.00   1.00  -102.09
#n.pH_p.tempZIP          6 221.14       26.44    0.00   1.00  -103.78
#n.pH_lcond_p.tempZIP    7 221.43       26.73    0.00   1.00  -102.64
#n.nat_lcond_p.tempZIP   7 221.99       27.30    0.00   1.00  -102.92
#n.cult_lcond_p.tempZIP  7 222.00       27.30    0.00   1.00  -102.92
#n.temp_pH_p.tempZIP     6 222.39       27.70    0.00   1.00  -104.40
#n.pH_OD_p.tempZIP       7 222.45       27.75    0.00   1.00  -103.15
#n.temp_p.tempZIP        5 222.77       28.07    0.00   1.00  -105.83
#NULO                    4 222.96       28.27    0.00   1.00  -107.12
#n.lives_pH_p.tempZIP    7 223.47       28.78    0.00   1.00  -103.66
#n.nat_pH_p.tempZIP      7 223.57       28.88    0.00   1.00  -103.71
#n.cult_pH_p.tempZIP     7 223.61       28.92    0.00   1.00  -103.73
#n.lives_p.tempZIP       6 223.77       29.08    0.00   1.00  -105.09
#n.OD_p.tempZIP          6 223.78       29.08    0.00   1.00  -105.10
#n.nat_p.tempZIP         6 223.94       29.24    0.00   1.00  -105.18
#n.cult_p.tempZIP        6 223.94       29.25    0.00   1.00  -105.18
#n.nat_temp_p.tempZIP    6 224.64       29.95    0.00   1.00  -105.53
#n.cult_temp_p.tempZIP   6 224.80       30.10    0.00   1.00  -105.61
#n.temp_OD_p.tempZIP     6 224.88       30.18    0.00   1.00  -105.65
#n.lives_temp_p.tempZIP  6 224.97       30.27    0.00   1.00  -105.69
#n.nat_OD_p.tempZIP      7 225.87       31.17    0.00   1.00  -104.86
#n.cult_OD_p.tempZIP     7 225.89       31.19    0.00   1.00  -104.87
#n.lives_OD_p.tempZIP    7 226.03       31.34    0.00   1.00  -104.94
#


###n.cult_elev_p.ZIP

#Call:
#pcount(formula = ~1 ~ elev + cult_500, data = ambyUMF, mixture = "ZIP")

#Abundance:
#  Estimate    SE    z  P(>|z|)
#(Intercept)    0.534 0.712 0.75 4.53e-01
#elev           2.229 0.340 6.56 5.45e-11
#cult_500       0.428 0.105 4.07 4.60e-05

#Detection:
#  Estimate    SE     z P(>|z|)
#-1.65 0.668 -2.46  0.0137

#Zero-inflation:
#  Estimate    SE       z P(>|z|)
#-0.0275 0.409 -0.0672   0.946

(confint(n.cult_elev_p.ZIP, type="state")) 

##############grfica de variables individualmente####################

png("psi.dry.count_col.png",
    width=16,height=16,units="cm",
    pointsize=10,bg="white",res=300)

par(mfrow = c(2,2), mar = c(5,5,2,2), cex.axis = 1.3, cex.lab = 1.5)


###################graficamos los valores predichos para cultivos#####

#se genera una nueva base
newcult <- data.frame(cult_500 = seq(min(scale(amby$cult_500)),
                                     max(scale(amby$cult_500)), 
                                     length=100),
                      elev = 0)

##con la cual se estiman los valores predichos
predcult <- predict(n.cult_elev_p.tempZIP, type = "state", newdata = newcult,
                    c.hat = 1.34)

#generamos una base de datos para graficar la respuesta a cultivos
cult.graph <- seq(min(amby$cult_500), max(amby$cult_500), length=100)

#extraemos los valores de los IC

cult <- (predcult$Predicted)

blower <- predcult$lower

bupper <- predcult$upper

###graficamos los valores predichos para cultivos a 500m modelo NB


plot(cult.graph, cult, type="l", las= 1, ylim = c(0,9.63),
     cex.axis = 1.5, cex.lab = 2, xlim = c(0,0.482),# family = "serif",
     xlab= "Crops 500m",
     ylab= "Expected abundance", bty="l", 
     col = "blue", lwd = 2)

matlines(cult.graph, bupper, lwd=2, col = "darkgray")

matlines(cult.graph, blower, lwd=2, col = "darkgray")


rug(amby$cult_500[amby$amby.t <="0"], side = 1, ticksize = 0.03, col = "black")

rug(amby$cult_500[amby$amby.t/amby$amby.t =="1"], side = 1, ticksize = 0.05, col = "blue")



#text(0.07, 9, "A)", family = "serif", cex = 3)

#####valores abundancia del modelo para la variable elevación

#se genera una nueva base para elevacion
newelev <- data.frame(elev = seq(min(scale(amby$elev)),
                                 max(scale(amby$elev)), 
                                 length=100),
                      cult_500 = 0)

##con la cual se estiman los valores predichos para modelo ZIP
predelev <- predict(n.cult_elev_p.tempZIP, type = "state",
                    newdata = newelev, c.hat =1.34 )

#generamos una base de datos para graficar la respuesta a lcond
elev.graph <- seq(min(amby$elev), max(amby$elev), length=100)

#extraemos los valores de los IC

elev <- (predelev$Predicted)

clower <- predelev$lower

cupper <- predelev$upper

###graficamos los valores predichos 

plot(elev.graph, elev, type="l", las= 1, ylim = c(0,144.4), xlim = c(500,2920),
     cex.axis = 1.5, cex.lab = 2,  #family = "serif",
     xlab= "Elevation (m asl)", ylab = "", bty="l",
     lwd = 2, col = "blue")

matlines(elev.graph, cupper, lwd = 2, col = "darkgray")

matlines(elev.graph, clower, lwd=2,  col = "darkgray")

rug(amby$elev[amby$amby.t <="0"], side = 1, ticksize = 0.03, col = "black")

rug(amby$elev[amby$amby.t/amby$amby.t =="1"], side = 1, ticksize = 0.05, col = "blue")



#text(700, 140, "B)", family = "serif", cex = 3)
#################covariables cultivo y elevación con el modelo binomial negativo##################
cl <- scale(seq(min(amby$cult_500), 
                max(amby$cult_500), length.out = 50)) # Standardised for prediction

ev <- scale(seq(min(amby$elev), 
                max(amby$elev), length.out = 50))


pred.matrix1 <- array(NA, dim = c(50, 50))#especificación de la matriz y sus dimenciones
for(i in 1:50){
  for(j in 1:50){
    newData <- data.frame(x=0, y=0, elev=ev[i], cult_500=cl[j], iLength=0)##ingreso de datos especificados para las predicciones conjuntas
    pred.matrix1[i,j] <- predict(n.cult_elev_p.tempZIP,type="state", 
                                 c.hat = 1.28, newdata=newData)[1,1]##generación de las predicciones a partir del modelo
  }
}


mapPalette <- colorRampPalette(c("gray","darkgrey", "yellow", "orange", "red"))

image(x=seq(0, 3000, length.out = 50), cex.lab = 2, cex.axis = 1.5,
      y=seq(0, 0.5, length.out = 50), 
      z=pred.matrix1, col = mapPalette(50), axes = F,
      xlab =  "Elevation (m asl)", 
      ylab =  "Crops 500m")

contour(x=seq(0, 3000, length.out = 50),
        y=seq(0, 0.5, length.out = 50),
        z=pred.matrix1, 
        add = T, col = "blue", labcex = 1.5,
        lwd = 2)

axis(1, at = seq(0, 3000, by = 500), cex.axis = 1.5)
axis(2, at = seq(0, 0.5, by = 0.1),  cex.axis = 1.5, las = 1)
box()

amby$amby.t <- rowSums(amby.y, na.rm = T)


points(amby$elev[amby$amby.t <= "0"], amby$cult_500[amby$amby.t <= "0"], pch = "+", cex = 1.5,
       col = "black")


linePallet <- colorRampPalette(c("white","orange","red"))

cl <- linePallet(17)

points(amby$elev[amby$amby.t > "0"], amby$cult_500[amby$amby.t > "0"], pch = "+", cex = 1.5,
       col = cl)


mtext("(a)", side = 3, line = -2, outer = T, at = 0.04, cex = 1.5,
      font = 1)

mtext("(b)", side = 3, line = -2, outer = T, at = 0.52, cex = 1.5,
      font = 1)


mtext("(c)", side = 3, line = -39, outer = T, at = 0.04, cex = 1.5,
      font = 1)

legend_image <- as.raster(matrix(mapPalette(5), nrow=1))

linePalletB <- colorRampPalette(c("black","white","orange","red"))

count_image <- as.raster(matrix(linePalletB(17), nrow=1))


plot(c(0,10),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x = 5, y = 0.95, labels = "Expected probability", cex= 1.5)
text(x= seq(2,8, l = 2), y =0.8, labels = c("Low", "High"), cex = 1.3)
rasterImage(legend_image, 2, 0.7, 8,0.6, angle = 0)

text(x = 5, y = 0.55, labels = "Estimated abundance", cex= 1.5)
text(x= seq(2,8, l = 2), y =0.44, labels = c("0", "17"), cex = 1.3)

rasterImage(count_image, 2, 0.4, 8,0.3, angle = 0)


legend(x = "bottom",
       legend=c("Estimated", "95% CI"),
       col=c("blue", "grey"),  lty=c(1,2), 
       cex=1.5, box.lty=0, text.font = 1, lwd = 2,
       xjust=0, yjust=0,text.width = 1.8, horiz = T,
       bg = NULL , x.intersp=0.3, seg.len = 0.7, y.intersp = 0.5)


dev.off()

save.image("abundancias.secas040719.RData")

amby.y2 <- amby[,2:7]

pseudoamby <- unmarkedFramePCount(y = amby.y2)

summary(ambyUMF)

plot(pseudoamby, col.regions = colorRampPalette(c('gray', 'darkred')))


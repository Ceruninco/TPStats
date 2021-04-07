library(maps)
library(TSPpackage)
library(sp)
library(Rcpp)
library(microbenchmark)
library(multcomp)
library(reshape)

villes <- read.csv('DonneesGPSvilles.csv',header=TRUE,dec='.',sep=';',quote="\"",fileEncoding="UTF-8")
str(villes)
coord <- cbind(villes$longitude,villes$latitude)
dist <- distanceGPS(coord)
voisins <- TSPnearest(dist)
print(voisins)
optCost <- TSPsolve(dist,'branch')
pathOpt <- c(1,8,9,4,21,13,7,10,3,17,16,20,6,19,15,18,11,5,22,14,12,2)

n <- 10

len <- 50
longueurs <- data.frame()
for (i in 1:len){
  sommets <- data.frame(x = runif(n), y = runif(n))
  couts <- distance(sommets)
  partial <- data.frame(x1 = TSPsolve(couts,'repetitive_nn'), x2 = TSPsolve(couts,'nearest_insertion'), x3 = TSPsolve(couts,'two_opt'), x4 = TSPsolve(couts,'nearest'), x5= TSPsolve(couts,'branch'))
  longueurs <- rbind(test,partial)
}
colnames(longueurs) <- c('repetitive_nn','nearest_insertion','two_opt','nearest','branch')

boxplot(longueurs)
#branch très bien

t.test(longueurs$repetitive_nn-longueurs$branch,alternative="greater")
#équivalent à
t.test(x=longueurs$repetitive_nn,y=longueurs$branch,paired=TRUE,alternative="greater")

data <- melt(longueurs) 

pairwise.t.test(data$value, data$variable, p.adjust.method='bonferroni', paired=TRUE)

microbenchmark(TSPsolve(x,'repetitive_nn'),TSPsolve(x,'nearest_insertion'), TSPsolve(x,'two_opt'), TSPsolve(x,'nearest'), TSPsolve(x,'branch'),
               times=50, setup={ 
                 sommets <- data.frame(x = runif(n), y = runif(n))
                 x <- distance(sommets)
               })

#partie 2
seqn <- seq(4,20,1)
times <- data.frame()
lenSeq <- length(seqn)

for(i in 1:lenSeq){
  partial <- microbenchmark(TSPsolve(couts, method = 'branch'),
                 times = 50,
                 setup = { n <- seqn[i]
                 couts <- distance(cbind(x = runif(n), y = runif(n)))}
  )$time
  times <- rbind(times,partial)
}

           
#semble exponentiel
par(mfrow=c(1,2)) # 2 graphiques sur 1 ligne
matplot(seqn, times, xlab='n', ylab='temps')
matplot(log(seqn), log(times), xlab='n', ylab=expression(log(times)))

vect_temps <- log(as.vector(times))
vect_temps <- melt(vect_temps)$value
vect_dim <- rep(log(seqn),times=50)
temps.lm <- lm(vect_temps~vect_dim)
summary(temps.lm)

par(mfrow=c(2,2)) # 4 graphiques, sur 2 lignes et 2 colonnes
plot(temps.lm)

shapiro.test(residuals(temps.lm))

#reg temps en fonction du nb de sommets
times.moy <- rowMeans(times)
vect_dim_2 <- log(seqn)^4
temps2.lm <- lm(times.moy~vect_dim_2)
summary(temps2.lm)
shapiro.test(residuals(temps2.lm))

#
data.graph <- read.csv('DonneesTSP.csv',header=TRUE,dec='.',sep=',',quote="\"",fileEncoding="UTF-8")
data.time <- log(unlist(subset(data.graph,select=c('tps'))))
data.features <- subset(data.graph,select=-c(tps))
data.features$dim <- log(data.features$dim)
data.lm <- lm(data.time~., data = data.features)
summary(data.lm)
shapiro.test(residuals(data.lm))

step(data.lm)

reducedModel <- lm(formula = data.time ~ dim + mean.long + sd.dist + mean.deg + 
     sd.deg + diameter, data = data.features)
shapiro.test(residuals(reducedModel))

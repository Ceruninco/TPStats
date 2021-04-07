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


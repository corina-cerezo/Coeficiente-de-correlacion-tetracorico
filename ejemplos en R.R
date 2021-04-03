library(psych)
?tetrachoric
tetrachoric(matrix(c(44268,193,14,0),2,2))  #MPLUS reports.24

tetrachoric(matrix(c(1562,42,383,94),2,2))  #Ejemplo de Ekrtöm
rt <- tetrachoric(matrix(c(631,125,147,147),2,2))  #ilustración 1 de Pearson
rt
rt$rho
rt$tau
chisq.test(matrix(c(1562,42,383,94),2,2))
qnorm((778-272)/(2*1050)+0.5)
qnorm(3/4)

rt <- tetrachoric(matrix(c(1766,842,842,722),2,2))  #ilustración 3 de Pearson
rt$rho
rt$tau
rt <- tetrachoric(matrix(c(254,136,156,193),2,2))  #ilustración 4 de Pearson
rt
rt$rho
rt$tau
rt$objective
rt$Call
data(bock)
responses <- table2df(bock.table[,2:6],count=bock.table[,7],
                      labs= paste("lsat6.",1:5,sep=""))
describe(responses)


# funcion para calcular el coeficiente de correlación tetracórico ---------
library(rootSolve)

rhot <- function(a,b,c,d){
  set.seed(2021)
  a = ifelse(a==0,0.5,a)
  b = ifelse(b==0,0.5,b)
  c = ifelse(c==0,0.5,c)
  d = ifelse(d==0,0.5,d)
  N = a+b+c+d
  h = qnorm(0.5 + ((a+c)-(b+d))/(2*N))
  k = qnorm(0.5 + ((a+b)-(c+d))/(2*N))
  H = (1/sqrt(2*pi))*exp(-0.5*h^2)
  K = (1/sqrt(2*pi))*exp(-0.5*k^2)
  izq = (a*d-b*c)/(N*N*H*K)
  coef1 = 1
  coef2 = h*k/factorial(2)
  coef3 = (h^2-1)*(k^2-1)/factorial(3)
  coef4 = h*(h^2-3)*k*(k^2-3)/factorial(4)
  coef5 = (h^4-6*h^2+3)*(k^4-6*k^2+3)/factorial(5)
  coef6 = h*(h^4-10*h^2+15)*k*(k^4-10*k^2+15)/factorial(6)
  coef7 = (h^6-15*h^4+45*h^2-15)*(k^6-15*k^4+45*k^2-15)/factorial(7)
  coef8 = h*(h^6-21*h^4+105*h^2-105)*k*(k^6-21*k^4+105*k^2-105)/factorial(8)
  serie <- function(x){
    return(-1*izq+coef1*x+coef2*x^2+coef3*x^3+coef4*x^4+coef5*x^5+coef6*x^6+coef7*x^7+coef8*x^8)
  }
  r <-  uniroot.all(serie, c(-1,1))
   
  beta1 <- (h-r*k)/sqrt(1-r^2)
  beta2 <- (k-r*h)/sqrt(1-r^2)
  psi1 <- pnorm(beta1)-0.5
  psi2 <- pnorm(beta2)-0.5
  chi0 <- (1/(2*pi))*(1/sqrt(1-r^2))*exp(-0.5*(1/(1-r^2))*(h^2+k^2-2*r*h*k))
  pe <- qnorm(3/4)/(sqrt(N)*chi0*N)*sqrt(((a+d)*(c+b))/(4)+psi2^2*((a+c)*(b+d))+psi1^2*((a+b)*(c+d))+2*psi1*psi2*(a*d-b*c)-psi2*(a*b-c*d)-psi1*(a*c-b*d))
  x <- data.frame(Estimacion = c(r),
                  P.E. = c(pe),
                  l.lim = c(r-pe),
                  u.lim = c(r+pe),
                  row.names = c("Coef. Corr."))
  
  return(x)
}
#Ilustracion 1 de Pearson
rhot(631,125,147,147)
#Ilustración 3 de Pearson
rhot(1766,842,842,722)
#Ilustración 6 de Pearson
rhot(1562,42,383,94)



# Error Probable ----------------------------------------------------------


# Mi tablita de frecuencias -----------------------------------------------
library(pracma)

frecuencias <- function(N,sigma1,sigma2,rhot,h,k){
  z <- function(x,y){
    N/(2*pi*sigma1*sigma2*sqrt(1-rhot^2))*exp(-0.5*(1/(1-rhot^2))*((x/sigma1)^2+(y/sigma2)^2-2*rhot*x*y/(sigma1*sigma2)))
  }
  a <- integral2(z, -10000, h, -10000, k)
  b <- integral2(z, h, 10000, -10000, k)
  c <- integral2(z, -10000, h, k, 10000)
  d <- integral2(z, h, 10000, k, 10000)
  return(c(round(a$Q),round(b$Q),round(c$Q),round(d$Q)))
}
#Intento 1
N <- 1000
s1 <- 5
s2 <- 5
p <- 0.3
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)

rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

library(psych)
library(pbivnorm)
library(rootSolve)



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
  peh <- qnorm(3/4)/(H*sqrt(N))*sqrt((b+d)*(a+c)/(N^2))
  pek <- qnorm(3/4)/(K*sqrt(N))*sqrt((c+d)*(a+b)/(N^2))
  eps = (a*d-b*c)/(N*N*H*K)
  coef1 = 1
  coef2 = h*k/factorial(2)
  coef3 = (h^2-1)*(k^2-1)/factorial(3)
  coef4 = h*(h^2-3)*k*(k^2-3)/factorial(4)
  coef5 = (h^4-6*h^2+3)*(k^4-6*k^2+3)/factorial(5)
  coef6 = h*(h^4-10*h^2+15)*k*(k^4-10*k^2+15)/factorial(6)
  coef7 = (h^6-15*h^4+45*h^2-15)*(k^6-15*k^4+45*k^2-15)/factorial(7)
  coef8 = h*(h^6-21*h^4+105*h^2-105)*k*(k^6-21*k^4+105*k^2-105)/factorial(8)
  serie1 <- function(x){
    return(-1*eps+coef1*x+coef2*x^2+coef3*x^3+coef4*x^4+coef5*x^5+coef6*x^6+coef7*x^7+coef8*x^8)
  }
  r1 <-  uniroot.all(serie1, c(-1,1))
  serie2 <- function(x){
    return(-1*eps+ x + h*k/2*x^2-(h^2+k^2-(h^2)*(k^2))/6*x^3 + h*k*((h^2)*(k^2)-3*(h^2+k^2)+5)/24*x^4)
  }
  tet <-  uniroot.all(serie2, c(-1,1))
  r2 <- sin(tet)
  #r2 <- ifelse(length(r2) == 0, "N.P.",r2)
  if(length(r1)==1 & length(r2)==1){
    per <- function (r){
      beta1 <- (h-r*k)/sqrt(1-r^2)
      beta2 <- (k-r*h)/sqrt(1-r^2)
      psi1 <- pnorm(beta1)-0.5
      psi2 <- pnorm(beta2)-0.5
      chi0 <- (1/(2*pi))*(1/sqrt(1-r^2))*exp(-0.5*(1/(1-r^2))*(h^2+k^2-2*r*h*k))
      return(qnorm(3/4)/(sqrt(N)*chi0*N)*sqrt(((a+d)*(c+b))/(4)+psi2^2*((a+c)*(b+d))+psi1^2*((a+b)*(c+d))+2*psi1*psi2*(a*d-b*c)-psi2*(a*b-c*d)-psi1*(a*c-b*d)))
    }
    x <- data.frame(Estimacion = c(r1,r2,h,k),
                    P.E. = c(per(r1),per(r2),peh,pek),
                    l.lim = c(r1-per(r1),r2-per(r2),h-peh,k-pek),
                    u.lim = c(r1+per(r1),r2+per(r2),h+peh,k+pek),
                    row.names = c("Coef. Corr. 1","Coef. Corr. 2","h","k"))
    
    return(x)
  }
  else{
    return(c("serie1:",r1,"la de sin:",r2))
  }
  
}
#Ilustración 1 de Pearson
rhot(631,125,147,147) #mi función
rt <- tetrachoric(matrix(c(631,125,147,147),2,2))  #función de la librería psych
rt
rt$rho
rt$tau

#Ilustración 3 de Pearson
rhot(1766,842,842,722)
#Ilustración 6 de Pearson
rhot(1562,42,383,94)

# Mi tablita de frecuencias -----------------------------------------------
library(pracma)
library(rmutil)

frecuencias <- function(N,sigma1,sigma2,rhot,h,k){
  z <- function(x,y){
    N/(2*pi*sigma1*sigma2*sqrt(1-rhot^2))*exp(-0.5*(1/(1-rhot^2))*((x/sigma1)^2+(y/sigma2)^2-2*rhot*x*y/(sigma1*sigma2)))
  }
  #Aquí ocupé una función para integrar pero falla con ciertos valores
  # a <- integral2(z, -10000,h, -10000, k)
  # b <- integral2(z, h, 10000, -10000, k)
  # c <- integral2(z, -10000, h, k, 10000)
  # d <- integral2(z, h, 10000, k, 10000)
  # return(c(round(a$Q,3),round(b$Q,3),round(c$Q,3),round(d$Q,3)))
  #Aquí ocupé otra función, corre bien pero es DEMASIADO LENTA
  # a <- int2(z, a=c(-10000,-10000), b=c(h,k))
  # b <- int2(z, a=c(h,-10000), b=c(10000,k))
  # c <- int2(z, a=c(-10000,k), b=c(h,10000))
  # d <- int2(z, a=c(h,k), b=c(10000,10000))
  #Aquí ocupé lo que debí haber hecho desde el principio
  aux1 <- pbivnorm::pbivnorm(x=c(h/s1), y=c(k/s2), rho=p)
  a <- N*aux1
  aux2 <- pbivnorm::pbivnorm(x=c(Inf), y=c(k/s2), rho=p)
  b <- N*(aux2-aux1)
  aux3 <- pbivnorm::pbivnorm(x=c(h/s1), y=c(Inf), rho=p)
  c <- N*(aux3-aux1)
  d <- N*(1-aux2-aux3+aux1)
  return(c(round(a),round(b),round(c),round(d)))
}
#Ejemplo 1
N <- 1000
s1 <- 5
s2 <- 5
p <- 0.4
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)
rhot(ns[1],ns[2],ns[3],ns[4])
M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

#Ejemplo 2, No se rechaza que p=0
N <- 10000
s1 <- 5
s2 <- 5
p <- 0
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)

rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

#Ejemplo 3, SI se rechaza que p=0 porque N es muy grande
N <- 1000000
s1 <- 5
s2 <- 5
p <- 0.05
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)

rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

#Ejemplo 4, NO se rechaza que p=0 porque N es pequeña
N <- 1000
s1 <- 5
s2 <- 5
p <- 0.05
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)

rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

#Ejemplo 5, OMG si corre si solo hay un cero (Como ya corregí frecuencias ya no hay ceros jo jo)
N <- 1000
s1 <- 5
s2 <- 5
p <- 0.49
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)
rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

rt <- tetrachoric(M)  #función de la librería psych
rt
rt$rho
rt$tau

#Ejemplo 6, como ya corregí las cosas no hay dos ceros xD
N <- 1000
s1 <- 5
s2 <- 5
p <- -0.89
h <- 4
k <- 1
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)

rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

rt <- tetrachoric(M)  #función de la librería psych
rt
rt$rho
rt$tau

#Ejemplo 7 ejemplo loco
N <- 1000
s1 <- 1
s2 <- 1
p <- -0.89
h <- 15 #vs 3
k <- 0.3
ns <- frecuencias(N,s1,s2,p,h,k)
ns
sum(ns)

rhot(ns[1],ns[2],ns[3],ns[4])

M <- as.table(rbind(c(ns[1],ns[2]), c(ns[3],ns[4])))
M
chisq.test(M)

rt <- tetrachoric(M)  #función de la librería psych
rt
rt$rho
rt$tau



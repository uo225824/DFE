library(orthopolynom)
library(nor1mix)

datos<- rnorMix(n, MW.nm3)
grado<-10
i<-10
aproximacion2<-function(datos,polinomio,criterio,hist,modelo){
  n1<-length(datos)
  grado<-n-1
  ll<-nmodes(datos,0.5,lowsup = -Inf,uppsup =Inf )
  Lt<-0
  for (i in 1:grado) {
    if (polinomio=="hermite") {
      
      a<-hermite.h.polynomials(i-1,normalized = T)
      phi<-as.matrix(as.data.frame(polynomial.values(polynomials=a,datos)))*exp(-datos^2/2)
      Zj<-as.matrix(as.data.frame(colMeans(phi)))
      S<-0
      for (i in 1:ncol(phi)) {
        S<-cbind(S,phi[,i]-Zj[i])
      }
      if (ncol(S)==2) {
        sigma2<-1/(n^2)*sum((S[,2:ncol(S)])^2)
        
      }else{
        sigma2<-1/(n^2)*colSums((S[,2:ncol(S)])^2)
      }
      sigma2<-as.matrix(sigma2)
      bj<-1-sigma2/Zj^2
      for (i in 1:length(bj)) {
        if (bj[i]<0) {
          bj[i]=0
        }
        if (bj[i]>1) {
          bj[i]=1
        }
      }
      the<-bj*Zj
      fr<-phi%*%the  
      for (j in 1:length(fr)) {
        if (fr[j]<0) {
          
            fr[j]=abs(fr[j])
          
        }
        #if (fr[j]<0) {
        #  fr=rep(0, length(fr))
        #}
      }
    }
    if (polinomio=="laguerre") {
      a<-laguerre.polynomials(i-1,normalized = T)
      phi<-as.matrix(as.data.frame(polynomial.values(polynomials=a,datos)))*exp(-datos/2)
      Zj<-as.matrix(as.data.frame(colMeans(phi)))
      S<-0
      for (i in 1:ncol(phi)) {
        S<-cbind(S,phi[,i]-Zj[i])
      }
      
      if (ncol(S)==2) {sigma2<-1/(n)*mean((S[,2:ncol(S)])^2)}
      else{sigma2<-1/(n)*colMeans((S[,2:ncol(S)])^2)}
      sigma2<-as.matrix(sigma2)
      bj<-(as.matrix(Zj^2)-sigma2)/Zj^2
      the<-bj*Zj
      fr<-phi%*%the 
      for (j in 1:length(fr)) {
        if (fr[j]<0) {
          fr=rep(0, length(fr))
        }
      }
      
    }
    if (criterio=="BIC") {
      L<--2*sum(log(fr))+i*log(n)
      Lt<-cbind(Lt,L)
    }
    if (criterio=="AICc") {
      if (kurtosis(datos)>0.5) {
        
        L<--2*sum(log(fr))+2*i+(2*i^2+2*i)/(n-i-1)-i*(kurtosis(datos))
        
      }else if (nmodes(datos,0.5,lowsup = -Inf,uppsup =Inf )>1) {
        
        L<--2*sum(log(fr))+2*i+(2*i^2+2*i)/(n-i-1)-ll*i}
      
      else{
        
        L<--2*sum(log(fr))+2*i+(2*i^2+2*i)/(n-i-1)
      }
    }

     
      Lt<-cbind(Lt,L)
    
    
  }
  Lt<-Lt[,2:ncol(Lt)]
  j1<-which.min(Lt);
  if (polinomio=="hermite") {
    a<-hermite.h.polynomials(j1,normalized = T)
    phi<-as.matrix(as.data.frame(polynomial.values(polynomials=a,datos)))*exp(-datos^2/2)
    y<-seq(min(datos),max(datos),0.001)
    Zj<-as.matrix(as.data.frame(colMeans(phi)))
    S<-0
    for (i in 1:ncol(phi)) {
      S<-cbind(S,phi[,i]-Zj[i])
    }
    
    if (ncol(S)==2) {sigma2<-1/(n^2)*sum((S[,2:ncol(S)])^2)}
    else{sigma2<-1/(n^2)*colSums((S[,2:ncol(S)])^2)}
    sigma2<-as.matrix(sigma2)
    bj<-1-sigma2/Zj^2
    the<-bj*Zj
    phi3<-as.matrix(as.data.frame(polynomial.values(polynomials=a,y)))*exp(-y^2/2)
    
    fr<-phi3%*%the
    datos<-sort(datos)
    if (hist=="si") {
      hist(datos,probability = T,ylim=c(0,1))
      lines(fr~y,col=3)
      lines(density(datos,bw="SJ"),col=2)
    }else{
      nms <- ls(pat="^MW.nm", "package:nor1mix")
      nms <- nms[order(as.numeric(substring(nms,6)))]
      nms<-nms[-17]
      plot(get(nms[modelo], "package:nor1mix"))
      
      lines(fr~y,col=3)
      lines(density(datos,bw="SJ"),col=2)}
    
    
  }
  if (polinomio=="laguerre") {
    a<-laguerre.polynomials(j1,normalized = T)
    phi<-as.matrix(as.data.frame(polynomial.values(polynomials=a,datos)))*exp(-datos/2)
    y<-seq(min(datos)-1,max(datos)+1,0.001)
    Zj<-as.matrix(as.data.frame(colMeans(phi)))
    S<-0
    for (i in 1:ncol(phi)) {
      S<-cbind(S,phi[,i]-Zj[i])
    }
    
    if (ncol(S)==2) {sigma2<-1/(n)*mean((S[,2:ncol(S)])^2)}
    else{sigma2<-1/(n)*colMeans((S[,2:ncol(S)])^2)}
    sigma2<-as.matrix(sigma2)
    bj<-(as.matrix(Zj^2)-sigma2)/Zj^2
    for (i in 1:length(bj)) {
      if (bj[i]>1) {
        bj[i]=1
      }
      if (bj[i]<0) {
        bj[i]=0
      }
    }
    the<-bj*Zj
    phi3<-as.matrix(as.data.frame(polynomial.values(polynomials=a,y)))*exp(-y/2)
    
    fr<-phi3%*%the
    datos<-sort(datos)
    
    if (hist=="si") {
      hist(datos,probability = T)
      lines(fr~y,col=3)
      lines(density(datos,bw="SJ"),col=2)
    }
    else{
      nms <- ls(pat="^MW.nm", "package:nor1mix")
      nms <- nms[order(as.numeric(substring(nms,6)))]
      nms<-nms[-17]
      plot(get(nms[modelo], "package:nor1mix"))
      
      lines(fr~y,col=3)
      lines(density(datos,bw="SJ"),col=2)}
    
    
  }
  return(cbind(j1,Lt))
}  

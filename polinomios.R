library(orthopolynom)
library(nor1mix)
library(multimode)

datos<- rnorMix(n, MW.nm3)
aproximacion<-function(datos,polinomio,criterio,hist,modelo){
  n<-length(datos)
  grado<-150
  ll<-nmodes(datos,0.5,lowsup = -Inf,uppsup =Inf )
  Lt<-0
  for (i in 1:grado) {
    if (polinomio=="hermite") {
      
    a<-hermite.h.polynomials(i-1,normalized = T)
    phi<-as.matrix(as.data.frame(polynomial.values(polynomials=a,datos)))*exp(-datos^2/2)
    medias<-as.matrix(as.data.frame(colMeans(phi)))
    fr<-phi%*%medias
    for (j in 1:length(fr)) {

      fr[j]=abs(fr[j])
    }
    }
    if (polinomio=="laguerre") {
      a<-laguerre.polynomials(i-1,normalized = T)
      phi<-as.matrix(as.data.frame(polynomial.values(polynomials=a,datos)))*exp(-datos/2)
      medias<-as.matrix(as.data.frame(colMeans(phi)))
      fr<-phi%*%medias
      for (j in 1:length(fr)) {
        if (fr[j]<0) {
          #fr=rep(0, length(fr))
          fr[j]=1
        }
      }
      
    }
    if (criterio=="BIC") {
      if (kurtosis(datos)>0.8) {
        
        L<--2*sum(log(fr))+i*log(n)-i*(kurtosis(datos))
        
      }else if (nmodes(datos,0.5,lowsup = -Inf,uppsup =Inf )>1) {
        
        L<--2*sum(log(fr))+i*log(n)-ll*i}
      
      else{
        
        L<--2*sum(log(fr))+i*log(n)}
      
      Lt<-cbind(Lt,L)
    }
    
    if (criterio=="AICc") {
      if (kurtosis(datos)>0.5) {
        
        L<--2*sum(log(fr))+2*i+(2*i^2+2*i)/(n-i-1)-i*(kurtosis(datos))
        
      }else if (nmodes(datos,0.5,lowsup = -Inf,uppsup =Inf )>3) {
        
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
    y<-seq(min(datos)-1,max(datos)+1,0.001)
    medias<-as.matrix(as.data.frame(colMeans(phi)))                
    phi3<-as.matrix(as.data.frame(polynomial.values(polynomials=a,y)))*exp(-y^2/2)
    
    fr<-phi3%*%medias
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
    if(abs(min(datos))>max(datos)){
      y<-seq(min(datos),abs(min(datos)),0.001)
    }else{
      y<-seq(-max(datos),max(datos),0.001)
    }
    
    medias<-as.matrix(as.data.frame(colMeans(phi)))                
    phi3<-as.matrix(as.data.frame(polynomial.values(polynomials=a,y)))*exp(-y/2)
    
    fr<-phi3%*%medias
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
  
#......................................................................#
#   Análise de diagnóstico do livro de Montegomery at al(outros autores)
#   Capítulo 7

#.........Concentração da madeira............#
x<-c(1.0, 1.5, 2.0, 3.0, 4.0, 4.5, 5, 
     5.5, 6.0, 6.5, 7.0, 8.0, 9:15)

#..........Tesnão do papel...................#
y<-c(06.3,11.1,20.0,24.0,26.1,30.0,33.8,34.0,38.1,39.9,
     42.0,46.1,53.1,52.0,52.5,48.0,42.8,27.8,21.9)


dados<-data.frame(x,y)
attach(dados)


#.....................................................#
require(ggplot2)

# Plotado o modleo poly de segunda ordem.
# Mas aparentemnte o de terceira ordem se adequa melhor
#.....................................................#


fit.model<-lm(y~I(x-mean(x))+ I((x-mean(x))^2), data = dados )
require(xtable)

xtable(summary(fit.model))

# Gráfico de dispersão com curva
ggplot(data = dados, aes(x,y))+
  geom_point()+
  geom_smooth(formula = y~poly(x,2), method = "lm")+
  labs(title = "Dispersão e polinômio de regressão estimado",
       x = "Hardwood Concentration (%)",
       y = "Tensile Strength (psi)")

# ou #

plot(x,y,pch=16, xlab = "Hardwood Concentration (%)",
     ylab = "Tensile Strength (psi)",
     main = "Dispersão e polinômio de regressão estimado") 

f<-function(x){fit.model$coefficients[1]+
    fit.model$coefficients[2]*(x-mean(x))+
    fit.model$coefficients[3]*(x-mean(x))^2}
curve(f(x), 0,16, col = "red",add = T)

#................Diagnóstico.................#
#............Resíduos studentizados..........#

X<-as.matrix(cbind(1,x,x^2))  # Matriz de especificação
H<-X%*%solve(t(X)%*%X)%*%t(X) 

n<-nrow(X)
p<-ncol(X)

SQRes<-t(fit.model$residuals)%*%fit.model$residuals
QMR<-SQRes/(n-p-1)

for(i in 1:n){
  e_studentizado<-fit.model$residuals/(sqrt(QMR*(1-H[i,i])))
}

#....Residuos studentizados vs valores ajustados  .....#
# Verificar se estão disposto aleatoriamente sobre 0 #

par(mfrow=c(2,1))
plot(fit.model$fitted.values,e_studentizado, pch = 16,
     xlab = "Valores ajustados", ylab = "Resíduo studentizados")
abline(0,0, col ="red")

#.........Resíduos studentizado externamentes.......#
#..Sem a presença do i-ésimo indivíduo discrepante..#

ti2<-e_studentizado*sqrt( (n-p-1)/(n - p - e_studentizado^2) )
plot(fit.model$fitted,ti2, 
     pch = 16,
     xlab = "Valores ajustados",
     ylab = "Resíduo studentizados Externamente")

abline(0,0, col ="red")
#===============================================================#
#....Análise da suposição de normalidade....#
# pag. 229 livro do Juvêncio at al
plot(fit.model)

diag.norm <- function(modelo=fit.model,iden=c(0,0,0,0,0,0),nome=seq(along = model.matrix(modelo)[,1])) {
  
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referências:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2ª ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4ª ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
  #
  
  if( class(modelo)[1]=="lm" || (class(modelo)[1]=="glm" && (modelo$family[[1]]=="Gaussian" | modelo$family[[1]]=="gaussian")) ) {
    
  } else {
    stop(paste("\nA classe do objeto deveria ser lm ou glm (com distribuicao gaussian) !!!"))
  }
  
  if(length(iden)<6) {
    iden<-c(-1,-1,-1,-1,-1,-1)
  }
  
  X <- model.matrix(modelo)
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  
  #para evitar divisão por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  lms <- summary(modelo)
  s <- lms$sigma
  r <- resid(modelo)
  ts <- r/(s*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ts^2)
  si <- lm.influence(modelo)$sigma
  tsi <- r/(si*sqrt(1-h))
  #dff <- sqrt(h/(1-h))*abs(tsi) #DFFIT
  ci <- sqrt( ((n-p)*h) / (p*(1-h)) )*abs(tsi) #aperfeiçoamento do DFFIT
  A <- diag(r)%*%H%*%diag(r)
  dmax <- abs(eigen(A)$vec[,1]/sqrt(eigen(A)$val[1]))
  m <- fitted(modelo)
  
  par(mfrow=c(2,3))
  
  plot(m,di,xlab="Valor Ajustado", ylab="Distância de Cook",main="Influência na Posição", ylim=c(0,max(di,2*mean(di))), pch=16)
  abline(2*mean(di),0,lty=2)
  while ( (!is.numeric(iden[1])) || (round(iden[1],0) != iden[1]) || (iden[1] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[1]<-as.numeric(out)
  }
  if(iden[1]>0) {identify(m,di,n=iden[1],labels=nome)}
  
  plot(m,ci,xlab="Valor Ajustado", ylab="Distância de Cook Modificada",main="Influência Posição/Escala", ylim=c(0,max(ci,2*mean(ci))), pch=16)
  abline(2*mean(ci),0,lty=2)
  while ( (!is.numeric(iden[2])) || (round(iden[2],0) != iden[2]) || (iden[2] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[2]<-as.numeric(out)
  }
  if(iden[2]>0) {identify(m,ci,n=iden[2],labels=nome)}
  
  plot(m,dmax,xlab="Valor Ajustado", ylab="dmax",main="Influência Local", ylim=c(0,max(dmax,2*mean(dmax))), pch=16)
  abline(2*mean(dmax),0,lty=2)
  while ( (!is.numeric(iden[3])) || (round(iden[3],0) != iden[3]) || (iden[3] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[3]<-as.numeric(out)
  }
  if(iden[3]>0) {identify(m,dmax,n=iden[3],labels=nome)}
  
  plot(m,h,xlab="Valor Ajustado", ylab="Medida h",main="Pontos Alavanca", ylim=c(0,max(h,2*p/n)), pch=16)
  abline(2*p/n,0,lty=2)
  while ( (!is.numeric(iden[4])) || (round(iden[4],0) != iden[4]) || (iden[4] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[4]<-as.numeric(out)
  }
  if(iden[4]>0) {identify(m,h,n=iden[4],labels=nome)}
  
  plot(m,tsi,xlab="Valor Ajustado", ylab="Resíduo Padronizado",main="Pontos Aberrantes", ylim=c(min(tsi)-1,max(tsi)+1), pch=16)
  abline(qt(.025,n-p-1),0,lty=2)
  abline(qt(1-.025,n-p-1),0,lty=2)
  while ( (!is.numeric(iden[5])) || (round(iden[5],0) != iden[5]) || (iden[5] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[5]<-as.numeric(out)
  }
  if(iden[5]>0) {identify(m,tsi,n=iden[5],labels=nome)}
  
  plot(m,abs(tsi),xlab="Valor Ajustado", ylab="Resíduo Padronizado Absoluto",main="Função de Variância", pch=16)
  lines(lowess(m,abs(tsi)))
  while ( (!is.numeric(iden[6])) || (round(iden[6],0) != iden[6]) || (iden[6] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[6]<-as.numeric(out)
  }
  if(iden[6]>0) {identify(m,abs(tsi),n=iden[6],labels=nome)}
  
  par(mfrow=c(1,1))
  list(ResPearsonStd=tsi,Di=di,Ci=ci,Dmax=dmax,h=h)
}

diag.norm(fit.model,iden=c(2,1,3,1,0,5))

#Clicar nos pontos para destacar
diag.norm(fit.model, iden = 2)

#................................................................#
#.................#......................#.......................#

envel.norm <- function(modelo=fit.model,iden=0,nome=seq(along = model.matrix(modelo)[,1]),sim=100,conf=.90,res=T,quad=T) {
  
  
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referências:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2ª ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4ª ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
  
  if( class(modelo)[1]=="lm" || (class(modelo)[1]=="glm" && (modelo$family[[1]]=="Gaussian" | modelo$family[[1]]=="gaussian")) ) {
    
  } else {
    stop(paste("\nA classe do objeto deveria ser lm ou glm (com distribuicao gaussian) !!!"))
  }
  
  alfa<-(1-conf)/2
  X <- model.matrix(modelo)
  y<-predict(modelo)+resid(modelo)
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  m <- fitted(modelo)
  
  #para evitar divisão por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  si <- lm.influence(modelo)$sigma
  r <- resid(modelo)
  tsi <- r/(si*sqrt(1-h))
  sigma<-summary(modelo)$sigma
  ti <- r/(sigma*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ti^2)
  
  e <- matrix(0,n,sim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:sim) {
    resp <- rnorm(n,m,sigma)
    fit <- lm(resp~X-1)
    ti<-resid(fit)/(summary(fit)$sigma*sqrt(1-h))
    if(res==F) {
      e[,i] <- (1/p)*(h/(1-h))*(ti^2)
    } else {
      e[,i] <- ti*sqrt( (n-p-1)/(n-p-(ti^2)) )
    }	
    e[,i] <- sort(e[,i])
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo,alfa)
    e2[i] <- quantile(eo,1-alfa)
  }
  
  med <- apply(e,1,median)
  
  if(quad==T) {
    par(pty="s")
  }
  if(res==F) {
    #Segundo McCullagh e Nelder (1989, pág.407) e Paula (2003, pág.57) deve-se usar qnorm((n+1:n+.5)/(2*n+1.125))
    #Segundo Neter et alli (1996, pág.597) deve-se usar qnorm((n+1:n-.125)/(2*n+0.5))
    qq<-qnorm((n+1:n+.5)/(2*n+1.125))
    plot(qq,sort(di),xlab="Quantil Meio-Normal",ylab="Distância de Cook", ylim=range(di,e1,e2), pch=16)
    nome<-nome[order(di)]
    r<-sort(di)
  } else {
    qq<-qnorm((1:n-.375)/(n+.25))
    plot(qq,sort(tsi),xlab="Quantil da Normal Padrão",ylab="Resíduo Studentizados", ylim=range(tsi,e1,e2), pch=16,main="Gráfico de envelope simulado com 90% de confiança")
    nome<-nome[order(tsi)]
    r<-sort(tsi)
  }
  lines(qq,e1,lty=1)
  lines(qq,e2,lty=1)
  lines(qq,med,lty=2)
  while ( (!is.numeric(iden)) || (round(iden,0) != iden) || (iden < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden<-as.numeric(out)
  }
  if(iden>0) {identify(qq,r,n=iden,labels=nome)}
  if(quad==T) {
    par(pty="m")
  }
  cat("\nBanda de ",conf*100,"% de confianca, obtida por ",sim," simulacoes.\n")
}

# Criando um gráfico de envelopes com 90% de confiança
envel.norm(fit.model,sim=10000,conf=.90)

# Identificando pontos
envel.norm(fit.model, iden = 4)

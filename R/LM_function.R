#' Analysis: Linear regression graph
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description Linear regression analysis of an experiment with a quantitative factor or isolated effect of a quantitative factor
#' @param trat Numerical vector with treatments (Declare as numeric)
#' @param resp Numerical vector containing the response of the experiment.
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param grau degree of the polynomial (1,2 or 3)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param legend.position legend position (\emph{default} is "top")
#' @param cardinal defines the value of y considered extreme (\emph{default} considers 0 germination)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @return The function returns the coefficients and respective p-values; statistical parameters such as AIC, BIC, R2, VIF; cardinal and optimal temperature and the graph using ggplot2 with the equation.
#' @keywords regression linear
#' @note if the maximum predicted value is equal to the maximum x, the curve does not have a maximum point within the studied range. If the minimum value is less than the lowest point studied, disregard the value.
#' @export
#' @examples
#' library(seedreg)
#' data("aristolochia")
#' attach(aristolochia)
#' LM_model(trat,resp, grau = 3)

LM_model=function(trat,
                  resp,
                  ylab="Germination (%)",
                  error="SE",
                  xlab=expression("Temperature ("^"o"*"C)"),
                  grau=NA,
                  theme=theme_classic(),
                  cardinal=0,
                  legend.position="top",
                  scale="none"){
  requireNamespace("ggplot2")
  requireNamespace("crayon")
  if(is.na(grau)==TRUE){grau=1}
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean)
  mod=c()
  mod1=c()
  mod2=c()
  mod05=c()

  modm=c()
  mod1m=c()
  mod2m=c()
  mod05m=c()

  text1=c()
  text2=c()
  text3=c()
  text05=c()

  mods=c()
  mod1s=c()
  mod2s=c()
  mod05s=c()

  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  dose=tapply(trat, trat, mean)
  moda=lm(resp~trat)
  mod1a=lm(resp~trat+I(trat^2))
  mod2a=lm(resp~trat+I(trat^2)+I(trat^3))
  mod05a=lm(resp~trat+I(trat^0.5))

  mods=summary(moda)$coefficients
  mod1s=summary(mod1a)$coefficients
  mod2s=summary(mod2a)$coefficients
  mod05s=summary(mod05a)$coefficients

  modm=lm(media~dose)
  mod1m=lm(media~dose+I(dose^2))
  mod2m=lm(media~dose+I(dose^2)+I(dose^3))
  mod05m=lm(media~dose+I(dose^0.5))

  if(grau=="1"){r2=round(summary(modm)$r.squared,2)}
  if(grau=="2"){r2=round(summary(mod1m)$r.squared,2)}
  if(grau=="3"){r2=round(summary(mod2m)$r.squared,2)}
  if(grau=="0.5"){r2=round(summary(mod05m)$r.squared,2)}
  if(grau=="1"){s1=s <- sprintf("~~~y == %e %s %e*x ~~~~~ italic(R^2) == %0.2f",
                                coef(moda)[1],
                                ifelse(coef(moda)[2] >= 0, "+", "-"),
                                abs(coef(moda)[2]),
                                as.numeric(as.character(r2)))}
  if(grau=="2"){s2=s <- sprintf("~~~y == %e %s %e * x %s %e * x^2 ~~~~~ italic(R^2) ==  %0.2f",
                                coef(mod1a)[1],
                                ifelse(coef(mod1a)[2] >= 0, "+", "-"),
                                abs(coef(mod1a)[2]),
                                ifelse(coef(mod1a)[3] >= 0, "+", "-"),
                                abs(coef(mod1a)[3]),
                                as.numeric(as.character(r2)))}
  if(grau=="3"){s3=s <- sprintf("~~~y == %e %s %e * x %s %e * x^2 %s %0.e * x^3 ~~~~~ italic(R^2) == %0.2f",
                                coef(mod2a)[1],
                                ifelse(coef(mod2a)[2] >= 0, "+", "-"),
                                abs(coef(mod2a)[2]),
                                ifelse(coef(mod2a)[3] >= 0, "+", "-"),
                                abs(coef(mod2a)[3]),
                                ifelse(coef(mod2a)[4] >= 0, "+", "-"),
                                abs(coef(mod2a)[4]),
                                as.numeric(as.character(r2)))}
  if(grau=="0.5"){s05=s <- sprintf("~~~y == %e %s %e * x %s %e * x^0.5 ~~~~~ italic(R^2) ==  %0.2f",
                                coef(mod05a)[1],
                                ifelse(coef(mod05a)[2] >= 0, "+", "-"),
                                abs(coef(mod05a)[2]),
                                ifelse(coef(mod05a)[3] >= 0, "+", "-"),
                                abs(coef(mod05a)[3]),
                                as.numeric(as.character(r2)))}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=as.numeric(as.character(names(media))),
                   resp=media,
                   desvio)
  grafico=ggplot(data1,aes(x=trat,y=resp))+
    geom_errorbar(aes(ymin=resp-desvio, ymax=resp+desvio),width=0.5)+
    geom_point(aes(fill=as.factor(rep(1,length(resp)))),na.rm=T,
               size=4.5,color="black",shape=21)+
    theme+ylab(ylab)+xlab(xlab)
  if(grau=="1"){grafico=grafico+geom_smooth(method = "lm",se=F, na.rm=T, formula = y~x,size=0.8,color="black")}
  if(grau=="2"){grafico=grafico+geom_smooth(method = "lm",se=F, na.rm=T, formula = y~x+I(x^2),size=0.8,color="black")}
  if(grau=="3"){grafico=grafico+geom_smooth(method = "lm",se=F, na.rm=T, formula = y~x+I(x^2)+I(x^3),size=0.8,color="black")}
  if(grau=="0.5"){grafico=grafico+geom_smooth(method = "lm",se=F, na.rm=T, formula = y~x+I(x^0.5),size=0.8,color="black")}
  if(grau=="1"){grafico=grafico+
    scale_fill_manual(values="gray",label=c(parse(text=s1)),name="")}
  if(grau=="2"){grafico=grafico+
      scale_fill_manual(values="gray",label=c(parse(text=s2)),name="")}
  if(grau=="3"){grafico=grafico+
      scale_fill_manual(values="gray",label=c(parse(text=s3)),name="")}
  if(grau=="0.5"){grafico=grafico+
    scale_fill_manual(values="gray",label=c(parse(text=s05)),name="")}

  grafico=grafico+
    theme(text = element_text(size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          axis.title = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.text=element_text(size=12),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)
  if(scale=="log"){grafico=grafico+scale_x_log10()}
  moda=lm(resp~trat)
  mod1a=lm(resp~trat+I(trat^2))
  mod2a=lm(resp~trat+I(trat^2)+I(trat^3))
  mod05a=lm(resp~trat+I(trat^0.5))

  if(grau=="1"){
  models=mods
  r2=summary(modm)$r.squared
  aic=AIC(moda)
  bic=BIC(moda)
  vif=NA
  temp1=seq(min(trat),max(trat),length.out=10000)
  result=predict(moda,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  fa=temp1[result<=cardinal & temp1>maximo]
  if(length(fa)>0){maxl=max(temp1[result<=cardinal & temp1>maximo])}else{maxl=NA}
  fb=temp1[result<=cardinal & temp1<maximo]
  if(length(fb)>0){minimo=max(temp1[result<=cardinal & temp1<maximo])}else{minimo=NA}}

  if(grau=="2"){
  models=mod1s
  r2=summary(mod1m)$r.squared
  aic=AIC(mod1a)
  bic=BIC(mod1a)
  vif=car::vif(mod1a)
  temp1=seq(min(trat),max(trat),length.out=10000)
  result=predict(mod1a,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  fa=temp1[result<=cardinal & temp1>maximo]
  if(length(fa)>0){maxl=max(temp1[result<=cardinal & temp1>maximo])}else{maxl=NA}
  fb=temp1[result<=cardinal & temp1<maximo]
  if(length(fb)>0){minimo=max(temp1[result<=cardinal & temp1<maximo])}else{minimo=NA}}


  if(grau=="3"){
  models=mod2s
  r2=summary(mod2m)$r.squared
  aic=AIC(mod2a)
  bic=BIC(mod2a)
  vif=car::vif(mod2a)
  temp1=seq(min(trat),max(trat),length.out=10000)
  result=predict(mod2a,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  fa=temp1[result<=cardinal & temp1>maximo]
  if(length(fa)>0){maxl=max(temp1[result<=cardinal & temp1>maximo])}else{maxl=NA}
  fb=temp1[result<=cardinal & temp1<maximo]
  if(length(fb)>0){minimo=max(temp1[result<=cardinal & temp1<maximo])}else{minimo=NA}}

  if(grau=="0.5"){
    models=mod05s
    r2=summary(mod05m)$r.squared
    aic=AIC(mod05a)
    bic=BIC(mod05a)
    vif=car::vif(mod05a)
    temp1=seq(min(trat),max(trat),length.out=10000)
    result=predict(mod05a,newdata = data.frame(trat=temp1),type="response")
    maximo=temp1[which.max(result)]
    respmax=result[which.max(result)]
    result=round(result,0)
    fa=temp1[result<=cardinal & temp1>maximo]
    if(length(fa)>0){maxl=max(temp1[result<=cardinal & temp1>maximo])}else{maxl=NA}
    fb=temp1[result<=cardinal & temp1<maximo]
    if(length(fb)>0){minimo=max(temp1[result<=cardinal & temp1<maximo])}else{minimo=NA}}
  cat("\n")
  graphs=data.frame("Parameter"=c("optimum temperature",
                                  "Maximum response",
                           "Predicted maximum value",
                           "Predicted minimum value",
                           "AIC","BIC","r-squared"),
             "values"=c(maximo,
                        respmax,
                        maxl,
                        minimo,
                        aic,
                        bic,
                        r2))
  graficos=list("Coefficients"=models,
                "test"=graphs,
                "VIF"=vif,
                grafico)
  print(graficos)
}

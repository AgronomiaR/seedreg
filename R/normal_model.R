#' Analysis: Normal model
#'
#' @param trat Numerical or complex vector with treatments
#' @param resp Numerical vector containing the response of the experiment.
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param legend.position legend position (\emph{default} is c(0.3,0.8))
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param cardinal defines the value of y considered extreme (\emph{default} considers 0 germination)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @return The function returns the coefficients and respective p-values; statistical parameters such as AIC, BIC, pseudo-R2; cardinal and optimal temperatures and the graph using ggplot2 with the equation.
#' @details The model function for the normal model is:
#' \deqn{f(x) = a \epsilon^{-\frac{(x-b)^2)}{c^2}}}
#' @note if the maximum predicted value is equal to the maximum x, the curve does not have a maximum point within the studied range. If the minimum value is less than the lowest point studied, disregard the value.
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @examples
#' library(seedreg)
#' data("aristolochia")
#' attach(aristolochia)
#' normal_model(trat,resp)

normal_model=function(trat,
                      resp,
                      ylab="Germination (%)",
                      xlab=expression("Temperature ("^"o"*"C)"),
                      theme=theme_classic(),
                      error="SE",
                      legend.position="top",
                      cardinal=0,
                      r2="all",
                      scale="none"){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  requireNamespace("crayon")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  mod=nls(resp ~ a*exp(-1/2*(trat-b)^2/c^2),
          start=c(a=100,b=mean(trat),c=5))
  coef=summary(mod)
  a=coef$coefficients[,1][1]
  b=coef$coefficients[,1][2]
  c=coef$coefficients[,1][3]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100

  equation=sprintf("~~~y==%0.3e*exp(frac(-(x %s %0.3e)^2, %0.3e)) ~~~~~ italic(R^2) == %0.2f",
                   a,
                   ifelse(b >= 0, "+", "-"),
                   abs(b),
                   2*c^2,
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  x=preditos$x
  y=preditos$y
  s=equation
  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  graph=ggplot(data,aes(x=xmean,y=ymean))
  if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                               width=0.5)}
  graph=graph+
    geom_point(aes(color="black"),size=4.5,shape=21,fill="gray")+
    theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,color="black"),size=0.8)+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=12),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  temp1=seq(min(trat),max(trat),length.out=10000)
  result=predict(mod,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  result=round(result,0)
  fa=temp1[result<=cardinal & temp1>maximo]
  if(length(fa)>0){maxl=max(temp1[result<=cardinal & temp1>maximo])}else{maxl=NA}
  fb=temp1[result<=cardinal & temp1<maximo]
  if(length(fb)>0){minimo=max(temp1[result<=cardinal & temp1<maximo])}else{minimo=NA}
  aic=AIC(mod)
  bic=BIC(mod)
  graphs=data.frame("Parameter"=c("optimum temperature",
                                  "Maximum response",
                                  "Predicted maximum value",
                                  "Predicted minimum value",
                                  "AIC","BIC","r-squared","RMSE"),
                    "values"=c(maximo,
                               respmax,
                               maxl,
                               minimo,
                               aic,bic,r2,rmse))
  graficos=list("Coefficients"=coef,
                "values"=graphs,
                graph)
  print(graficos)

}

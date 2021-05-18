#' Analysis: graph for not significant trend
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description Ggraph for non-significant trend. Can be used within the multicurve command
#' @param trat Numerical vector with treatments (Declare as numeric)
#' @param resp Numerical vector containing the response of the experiment.
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param legend.position legend position (\emph{default} is "top")
#' @return The function returns an exploratory graph of segments
#' @keywords non-significant
#' @export
#' @examples
#' library(seedreg)
#' data("aristolochia")
#' attach(aristolochia)
#' N_model(trat,resp)

N_model=function(trat,
                  resp,
                  ylab="Germination (%)",
                  error="SE",
                  xlab=expression("Temperature ("^"o"*"C)"),
                  theme=theme_classic(),
                  legend.position="top"){
  requireNamespace("ggplot2")
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean)
  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=as.numeric(as.character(names(media))),
                   resp=media,
                   desvio)
  temp1=dose
  result=media
  s="not~significant"
  grafico=ggplot(data1,aes(x=trat,y=resp))+
    geom_errorbar(aes(ymin=resp-desvio, ymax=resp+desvio),width=0.5)+
    geom_point(aes(color=as.factor(rep(1,length(resp)))),na.rm=T,
               size=4.5,fill="gray",shape=21)+
    theme+ylab(ylab)+xlab(xlab)+
    scale_color_manual(values="gray",label=c(parse(text=s)),name="")+
    theme(text = element_text(size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          axis.title = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.text=element_text(size=12),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)
  graficos=list(teste="not significant",grafico)
  print(graficos)
}

#' Merge multiple curves into a single graph
#' @param plots list with objects of type LM_model, BC_model, CD_model, LL_model or normal_model
#' @param theme ggplot2 theme (\emph{default} is theme_classi())
#' @param legend.title caption title
#' @param trat name of the curves
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param legend.position legend position (\emph{default} is c(0.3,0.8))
#' @param gray gray scale (\emph{default} is FALSE)
#' @return The function returns a graph joining the outputs of the functions LM_model, LL_model, BC_model, CD_model, loess_model, normal_model, piecewise_model and N_model
#' @author Gabriel Danilo Shimizu
#' @export
#' @examples
#' library(seedreg)
#' data("aristolochia")
#' attach(aristolochia)
#' a=LM_model(trat,resp)
#' b=LL_model(trat,resp,npar = "LL.3")
#' c=BC_model(trat,resp, npar = "BC.4")
#' d=CD_model(trat,resp, npar = "CRS.4")
#' multicurve(list(a,b,c,d))

multicurve=function(plots,
               theme=theme_classic(),
               legend.title=NULL,
               legend.position="top",
               trat=NA,
               gray=FALSE,
               ylab="Germination (%)",
               xlab=expression("Temperature ("^"o"*"C)")){
  requireNamespace("ggplot2")
  equation=1:length(plots)
  grafico=ggplot()
  if(gray==FALSE){
  for(i in 1:length(plots)){
    equation[[i]]=plots[[i]][[]]$plot$s
    x=plots[[i]][[]]$plot$temp1
    y=plots[[i]][[]]$plot$result
    data=data.frame(x,y,color=as.factor(i))
    pontosx=plots[[i]][[]]$plot$data1$trat
    pontosy=plots[[i]][[]]$plot$data1$resp
    desvio=plots[[i]][[]]$plot$desvio
    pontos=data.frame(x=pontosx,
                      y=pontosy,
                      desvio=desvio,
                      color=as.factor(i))
    color=pontos$color
    grafico=grafico+
      geom_errorbar(data=pontos,
                    aes(x=x,
                        y=y,
                        ymin=y-desvio,
                        ymax=y+desvio,
                        color=color,
                        group=color),width=0.1, size=0.8)+
      geom_point(data=pontos,aes(x=x,y=y,
                                 color=color,
                                 group=color),size=4.5)+
      geom_line(data=data,aes(x=x,
                              y=y,
                              color=color,
                              group=color),size=0.8)
  }
  texto=parse(text=paste(trat,"~",unlist(equation)))
  grafico=grafico+
    scale_color_discrete(label=texto)+
    theme+labs(color=legend.title)+
    theme(axis.text = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.justification='left',
          legend.direction = "vertical",
          legend.text.align = 0)+ylab(ylab)+xlab(xlab)}
  if(gray==TRUE){
    for(i in 1:length(plots)){
      equation[[i]]=plots[[i]][[]]$plot$s
      x=plots[[i]][[]]$plot$temp1
      y=plots[[i]][[]]$plot$result
      data=data.frame(x,y,color=as.factor(i))
      pontosx=plots[[i]][[]]$plot$data1$trat
      pontosy=plots[[i]][[]]$plot$data1$resp
      desvio=plots[[i]][[]]$plot$desvio
      pontos=data.frame(x=pontosx,y=pontosy,desvio=desvio,color=as.factor(i))
      grafico=grafico+
        geom_errorbar(data=pontos,
                      aes(x=x,
                          y=y,
                          ymin=y-desvio,
                          ymax=y+desvio),width=0.1, size=0.8)+
        geom_point(data=pontos,aes(x=x,
                                   y=y,
                                   pch=color,
                                   group=color),
                   size=4.5,fill="gray")+
        geom_line(data=data,aes(x=x,
                                y=y,
                                lty=color,
                                group=color),size=0.8)
    }
    texto=parse(text=paste(trat,"~",unlist(equation)))
    grafico=grafico+
      scale_linetype_discrete(label=texto)+
      scale_shape_discrete(label=texto)+
      theme+labs(lty=legend.title,shape=legend.title)+
      theme(axis.text = element_text(size=12,color="black"),
            legend.position = legend.position,
            legend.text=element_text(size=12),
            legend.justification='left',
            legend.direction = "vertical",
            legend.text.align = 0)+ylab(ylab)+xlab(xlab)}
  print(grafico)
}

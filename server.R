source("functions.R")
server <- function(input, output) {

###########################################################################################################################################################################
### graph of distribution turning points
  
###########################################################################################################################################################################
  output$distributionplot <- renderPlot({
    distr <- as.numeric(unlist(strsplit(input$distr,",")))
    D=input$D
    barplot(distr~seq(0,D),xlab="Turning point",ylab="Proportion",ylim=c(0,1))
    })

###########################################################################################################################################################################
### dropout graphs
###########################################################################################################################################################################
  output$survivalplots <- renderPlot({
    omega=input$omega
    gamma=input$gamma
    D=input$D
    time=seq(0,1,length=D+1)
    remain = (1 - omega)^time^gamma
    par(mfrow=c(1,2))
    plot(seq(0,D),remain,type="b",pch=16,ylim=c(0,1),xlim=c(0,D),xlab="time",ylab="probability",main="Survival")
    hazard=rep(0,D)
    for(kk in 1:D)    
      hazard[kk]=(remain[kk]-remain[kk+1])/remain[kk]
    plot(seq(1,D),hazard,type="b",pch=16,xlim=c(0,D),xlab="time",ylab="probability", main="Hazard")
  })
  
###########################################################################################################################################################################
### Results graph 1: power versus sample size for mean growth rate 1
###########################################################################################################################################################################
  output$Resultsplot1 <- renderPlotly({
    
    validate(need(input$D>=5,"Duration of the study should be at least five"))
    validate(need(input$D==round(input$D),"Duration of the study should be an integer"))
    
    validate(need(input$var.e>=0,"Residual variance should be positive"))
    validate(need(input$var.u0>=0,"Variance random intercept should be positive"))
    validate(need(input$var.u1>=0,"Variance random slope 1 should be positive"))
    validate(need(input$var.u2>=0,"Variance random slope 2 should be positive"))
    
    validate(need(input$alpha>=0&input$alpha<=1,"Type I error rate should be between 0 and 1"))
    
    validate(need(input$omega>=0&input$omega<=1,"Parameter omega should be between 0 and 1"))
    validate(need(input$gamma>0,"Parameter gamma should be above 0"))
    
    distr <- as.numeric(unlist(strsplit(input$distr,",")))
    D=input$D
    validate(need(sum(distr)==1,"Sum of proportions in distribution should be equal to 1"))
    validate(need(length(distr)==D+1,"Length of distribution vector should be equal to D+1"))
    validate(need(distr[1]<1,"In the case all subjects have the same turning point: this turning point cannot occur at the first measurement occasion"))
    validate(need(distr[D+1]<1,"In the case all subjects have the same turning point: this turning point cannot occur at the last measurement occasion"))
    
    
    var.e=input$var.e
    var.u0=input$var.u0
    var.u1=input$var.u1
    var.u2=input$var.u2
    cov.u01=input$covar.u01
    cov.u02=input$covar.u02
    cov.u12=input$covar.u12
    DD=matrix(0,nrow=3,ncol=3)
    diag(DD)=c(var.u0,var.u1,var.u2)
    DD[1,2]=DD[2,1]=cov.u01
    DD[1,3]=DD[3,1]=cov.u02
    DD[2,3]=DD[3,2]=cov.u12
    
    N=seq(1,10000) # sample size
    m=input$D+1 # number of measurements (duration + 1 since also a measurement at baseline)
    omega=input$omega 
    gamma=input$gamma

    infmat.no=matrix(0,nrow=3,ncol=3) # information matrix if attrition absent
    for(ii in 0:D)
        {
        infmatrix=f.infmat(var.e,DD,N=1,m=m,T=ii,omega=0,gamma=gamma)
        infmat.no=infmat.no+distr[ii+1]*infmatrix
    }
    covmat.no=solve(infmat.no)
    varbeta.no=covmat.no[2,2]
      
    infmat.yes=matrix(0,nrow=3,ncol=3) # information matrix if attrition present
    for(ii in 0:D)
      {
      infmatrix=f.infmat(var.e,DD,N=1,m=m,T=ii,omega=omega,gamma=gamma)
      infmat.yes=infmat.yes+distr[ii+1]*infmatrix
    }
    covmat.yes=solve(infmat.yes)
    varbeta.yes=covmat.yes[2,2]

    beta1=abs(input$beta1)
    alpha=input$alpha
    test=input$test
    
    varbeta.no=varbeta.no/N
    if(test==1)
      power.no=pnorm(beta1/sqrt(varbeta.no)-qnorm(1-alpha))
    if(test==2)
      power.no=pnorm(beta1/sqrt(varbeta.no)-qnorm(1-alpha/2))
    
    varbeta.yes=varbeta.yes/N
    if(test==1)
      power.yes=pnorm(beta1/sqrt(varbeta.yes)-qnorm(1-alpha))
    if(test==2)
      power.yes=pnorm(beta1/sqrt(varbeta.yes)-qnorm(1-alpha/2))

    plotdata=cbind(N,power.no,power.yes)
    plotdata=plotdata[plotdata[,3]<0.99,]
    plotdata=as.data.frame(plotdata)

    key <- row.names(plotdata)
    plot_ly(x = ~plotdata[,1], y = ~plotdata[,2], name = 'No attrition', key = ~key, mode ='lines',type="scatter",width=1,showlegend=TRUE,line = list(color = 'green')) %>%
      add_trace(y = ~plotdata[,3], name = 'Attrition', key = ~key, mode = 'lines',line = list(color = 'red')) %>%
      layout(legend = list(orientation = "h", xanchor = "center", x = 0.5,y=-0.25)) %>% 
      layout(dragmode = "select",
             xaxis = list(title = "Sample size",  showgrid = F,zeroline=TRUE,showline = TRUE),
             yaxis = list(title = "Power", range=c(0,1),showgrid = T,zeroline=TRUE,showline = TRUE))
  })

###########################################################################################################################################################################
### Results graph 2: power versus sample size for mean growth rate 2
###########################################################################################################################################################################
  output$Resultsplot2 <- renderPlotly({
    
    validate(need(input$D>=5,"Duration of the study should be at least five"))
    validate(need(input$D==round(input$D),"Duration of the study should be an integer"))
    
    validate(need(input$var.e>=0,"Residual variance should be positive"))
    validate(need(input$var.u0>=0,"Variance random intercept should be positive"))
    validate(need(input$var.u1>=0,"Variance random slope 1 should be positive"))
    validate(need(input$var.u2>=0,"Variance random slope 2 should be positive"))
    
    validate(need(input$alpha>=0&input$alpha<=1,"Type I error rate should be between 0 and 1"))
    
    validate(need(input$omega>=0&input$omega<=1,"Parameter omega should be between 0 and 1"))
    validate(need(input$gamma>0,"Parameter gamma should be above 0"))
    
    distr <- as.numeric(unlist(strsplit(input$distr,",")))
    D=input$D
    validate(need(sum(distr)==1,"Sum of proportions in distribution should be equal to 1"))
    validate(need(length(distr)==D+1,"Length of distribution vector should be equal to D+1"))
    validate(need(distr[1]<1,"In the case all subjects have the same turning point: this turning point cannot occur at the first measurement occasion"))
    validate(need(distr[D+1]<1,"In the case all subjects have the same turning point: this turning point cannot occur at the last measurement occasion"))
    

    var.e=input$var.e
    var.u0=input$var.u0
    var.u1=input$var.u1
    var.u2=input$var.u2
    cov.u01=input$covar.u01
    cov.u02=input$covar.u02
    cov.u12=input$covar.u12
    DD=matrix(0,nrow=3,ncol=3)
    diag(DD)=c(var.u0,var.u1,var.u2)
    DD[1,2]=DD[2,1]=cov.u01
    DD[1,3]=DD[3,1]=cov.u02
    DD[2,3]=DD[3,2]=cov.u12
    
    N=seq(1,10000)
    m=input$D+1
    omega=input$omega
    gamma=input$gamma

    infmat.no=matrix(0,nrow=3,ncol=3) # information matrix if attrition absent
    for(ii in 0:D)
    {
      infmatrix=f.infmat(var.e,DD,N=1,m=m,T=ii,omega=0,gamma=gamma)
      infmat.no=infmat.no+distr[ii+1]*infmatrix
    }
    covmat.no=solve(infmat.no)
    varbeta.no=covmat.no[3,3]
    
    infmat.yes=matrix(0,nrow=3,ncol=3) # information matrix if attrition present
    for(ii in 0:D)
    {
      infmatrix=f.infmat(var.e,DD,N=1,m=m,T=ii,omega=omega,gamma=gamma)
      infmat.yes=infmat.yes+distr[ii+1]*infmatrix
    }
    covmat.yes=solve(infmat.yes)
    varbeta.yes=covmat.yes[3,3]
    
    beta1=abs(input$beta2)
    alpha=input$alpha
    test=input$test
    
    varbeta.no=varbeta.no/N
    if(test==1)
      power.no=pnorm(beta1/sqrt(varbeta.no)-qnorm(1-alpha))
    if(test==2)
      power.no=pnorm(beta1/sqrt(varbeta.no)-qnorm(1-alpha/2))
    
    varbeta.yes=varbeta.yes/N
    if(test==1)
      power.yes=pnorm(beta1/sqrt(varbeta.yes)-qnorm(1-alpha))
    if(test==2)
      power.yes=pnorm(beta1/sqrt(varbeta.yes)-qnorm(1-alpha/2))
    
    plotdata=cbind(N,power.no,power.yes)
    plotdata=plotdata[plotdata[,3]<0.99,]
    plotdata=as.data.frame(plotdata)

    key <- row.names(plotdata)
    plot_ly(x = ~plotdata[,1], y = ~plotdata[,2], name = 'No attrition', key = ~key, mode ='lines',type="scatter",width=1,showlegend=TRUE,line = list(color = 'green')) %>%
      add_trace(y = ~plotdata[,3], name = 'Attrition', key = ~key, mode = 'lines',line = list(color = 'red')) %>%
      layout(legend = list(orientation = "h", xanchor = "center", x = 0.5,y=-0.25)) %>% 
      layout(dragmode = "select",
             xaxis = list(title = "Sample size",  showgrid = F,zeroline=TRUE,showline = TRUE),
             yaxis = list(title = "Power", range=c(0,1),showgrid = T,zeroline=TRUE,showline = TRUE))
  })
  
###########################################################################################################################################################################
### Results graph 3: power versus sample size for difference in mean growth rates
###########################################################################################################################################################################
  output$Resultsplot3 <- renderPlotly({
    
    validate(need(input$D>=5,"Duration of the study should be at least five"))
    validate(need(input$D==round(input$D),"Duration of the study should be an integer"))
    
    validate(need(input$var.e>=0,"Residual variance should be positive"))
    validate(need(input$var.u0>=0,"Variance random intercept should be positive"))
    validate(need(input$var.u1>=0,"Variance random slope 1 should be positive"))
    validate(need(input$var.u2>=0,"Variance random slope 2 should be positive"))
    
    validate(need(input$alpha>=0&input$alpha<=1,"Type I error rate should be between 0 and 1"))
    
    validate(need(input$omega>=0&input$omega<=1,"Parameter omega should be between 0 and 1"))
    validate(need(input$gamma>0,"Parameter gamma should be above 0"))
    
    distr <- as.numeric(unlist(strsplit(input$distr,",")))
    D=input$D
    validate(need(sum(distr)==1,"Sum of proportions in distribution should be equal to 1"))
    validate(need(length(distr)==D+1,"Length of distribution vector should be equal to D+1"))
    validate(need(distr[1]<1,"In the case all subjects have the same turning point: this turning point cannot occur at the first measurement occasion"))
    validate(need(distr[D+1]<1,"In the case all subjects have the same turning point: this turning point cannot occur at the last measurement occasion"))
    

    var.e=input$var.e
    var.u0=input$var.u0
    var.u1=input$var.u1
    var.u2=input$var.u2
    cov.u01=input$covar.u01
    cov.u02=input$covar.u02
    cov.u12=input$covar.u12
    DD=matrix(0,nrow=3,ncol=3)
    diag(DD)=c(var.u0,var.u1,var.u2)
    DD[1,2]=DD[2,1]=cov.u01
    DD[1,3]=DD[3,1]=cov.u02
    DD[2,3]=DD[3,2]=cov.u12
    
    N=seq(1,10000)
    m=input$D+1
    omega=input$omega
    gamma=input$gamma

    contrast=matrix(c(0,1,-1),nrow=3,ncol=1)
    
    infmat.no=matrix(0,nrow=3,ncol=3)
    for(ii in 0:D)
    {
      infmatrix=f.infmat(var.e,DD,N=1,m=m,T=ii,omega=0,gamma=gamma)
      infmat.no=infmat.no+distr[ii+1]*infmatrix
    }
    covmat.no=solve(infmat.no)
    varcontrast.no=t(contrast)%*%covmat.no%*%contrast
    
    infmat.yes=matrix(0,nrow=3,ncol=3)
    for(ii in 0:D)
    {
      infmatrix=f.infmat(var.e,DD,N=1,m=m,T=ii,omega=omega,gamma=gamma)
      infmat.yes=infmat.yes+distr[ii+1]*infmatrix
    }
    covmat.yes=solve(infmat.yes)
    varcontrast.yes=t(contrast)%*%covmat.yes%*%contrast
    
    beta1=input$beta1
    beta2=input$beta2
    beta.diff=abs(beta1-beta2)
    alpha=input$alpha
    test=input$test
    
    varcontrast.no=as.vector(varcontrast.no)/N
    if(test==1)
      power.no=pnorm(beta.diff/sqrt(varcontrast.no)-qnorm(1-alpha))
    if(test==2)
      power.no=pnorm(beta.diff/sqrt(varcontrast.no)-qnorm(1-alpha/2))
    
    varcontrast.yes=as.vector(varcontrast.yes)/N
    if(test==1)
      power.yes=pnorm(beta.diff/sqrt(varcontrast.yes)-qnorm(1-alpha))
    if(test==2)
      power.yes=pnorm(beta.diff/sqrt(varcontrast.yes)-qnorm(1-alpha/2))

    plotdata=cbind(N,power.no,power.yes)
    plotdata=plotdata[plotdata[,3]<0.99,]
    plotdata=as.data.frame(plotdata)

    key <- row.names(plotdata)
    plot_ly(x = ~plotdata[,1], y = ~plotdata[,2], name = 'No attrition', key = ~key, mode ='lines',type="scatter",width=1,showlegend=TRUE,line = list(color = 'green')) %>%
      add_trace(y = ~plotdata[,3], name = 'Attrition', key = ~key, mode = 'lines',line = list(color = 'red')) %>%
      layout(legend = list(orientation = "h", xanchor = "center", x = 0.5,y=-0.25)) %>% 
      layout(dragmode = "select",
             xaxis = list(title = "Sample size",  showgrid = F,zeroline=TRUE,showline = TRUE),
             yaxis = list(title = "Power", range=c(0,1),showgrid = T,zeroline=TRUE,showline = TRUE))
  })

}
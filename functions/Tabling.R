# ======================================
#      Project: Tabling functions
#      Casey Crisman-Cox
#      Useful Functions
#      Tabling.R
# ======================================
#######Dependencies#####
if("xtable" %in% installed.packages()){
     library(xtable)
}else{
     install.packages("xtable")
     library(xtable)
}
if("stringr" %in% installed.packages()){
     library(stringr)
}else{
     install.packages("stringr")
     library(stringr)
}


num2str <- function(x){formatC(x, digits=2, format='f')}

mod.sum<-function(mod.list, 
                  apsr=FALSE, 
                  stars=c("default", "all", "CI"),
                  seNote="Standard Errors in Parenthesis",
                  model.numbers=NULL,
                  model.names=NULL,
                  dep.varnames=NULL,
                  k=1){  ##Function to add in LogLik and N to xtable
     if(is.null(model.numbers)){
          model.numbers <- paste("\\multicolumn{1}{c}{Model ",
                               k:(k+length(mod.list)),
                               "}",
                               sep="")
     }
     if(is.null(dep.varnames)){
          dep.varnames <- rep("%", length(mod.list))
     }
     if(is.null(model.names)){
          model.names <- rep("%", length(mod.list))
     }
     
     out<-NULL
     if(apsr==FALSE){
          mod<-mod.list[[1]]
          ##Works on maxLik, lm, glm, game, and surv
          if("maxLik" %in% class(mod)){ ##Check if the class is maxLik
               out<-paste("\\hline ", 
                          "Log $L$ & \\multicolumn{4}{c}{",num2str(logLik(mod)), "}\\\\\n",
                          "$N$ & \\multicolumn{4}{c}{", nrow(mod$gradientObs), "}\\\\\n", 
                          sep="") ##Creates two rows which we'll add to xtable
          }
          if("lm" == class(mod)[1]){
               out<-paste("\\hline ",
                          "adj. R$^2$& \\multicolumn{4}{c}{", num2str(summary(mod)$adj.r.squared),
                          "}\\\\\n",
                          "$N$ & \\multicolumn{4}{c}{", nrow(mod$model), "}\\\\\n",
                          sep="") ##Same as above
          }
          if("plm" == class(mod)[1]){
               out<-paste("\\hline ",
                          "adj. R$^2$& \\multicolumn{4}{c}{", summary(baseFETE)$r.squared[2],
                          "}\\\\\n",
                          "$N$ & \\multicolumn{4}{c}{", nrow(mod$model), "}\\\\\n",
                          sep="") ##Same as above
          }
          
          if("glm" %in% class(mod) | "game" %in% class(mod) 
             |"polr" %in% class(mod)){
               out<-paste("\\hline ",
                          "Log $L$ & \\multicolumn{4}{c}{",num2str(logLik(mod)), "}\\\\\n",
                          "$N$ & \\multicolumn{4}{c}{", nrow(mod$model), "}\\\\\n", 
                          sep="") ##Creates two rows which we'll add to xtable
          }
          if("polywog" == class(mod)[1]){
               out<-paste("\\hline ",
                          "Lambda & \\multicolumn{4}{c}{",num2str(mod$lambda), "}\\\\\n",
                          "$N$ & \\multicolumn{4}{c}{", mod$nobs, "}\\\\\n", 
                          sep="") ##Creates two rows which we'll add to xtable
          }
          if("list" == class(mod)[1]){
               if(is.null(mod$ar2)){
                    out<-paste("\\hline ",
                          "Log $L$ & \\multicolumn{4}{c}{",num2str(mod$logLik), "}\\\\\n",
                          "$N$ & \\multicolumn{4}{c}{", num2str(mod$Obs), "}\\\\\n", 
                          sep="") ##Creates two rows which we'll add to xtable
               }else{
                    out<-paste("\\hline ",
                               "adj. R$^2$& \\multicolumn{4}{c}{", num2str(mod$ar2),
                               "}\\\\\n",
                               "$N$ & \\multicolumn{4}{c}{", nrow(mod$model), "}\\\\\n",
                               sep="") ##Same as above
               }
          }
          out<-list(out=out, 
                    length=length(unique(names(coef(mod)))))
     }else{
          out.L<-c("Log $L$")
          out.A <- c("$\\alpha$")
          out.Lam<-c("$\\lambda$")
          out.Rho<-c("$\\rho$")
          out.R<-c("adj. R$^2$")
          out.N<-c("$N$")
          modNums <- " "
          if(all(dep.varnames=="%")){ 
               depNames <- "%"
          }else{
               depNames <- " "
          }
          if(all(model.names=="%")){ 
               modNames <- "%"
          }else{
               modNames <- " "
          }
          
          for(i in 1:length(mod.list)){
               mod.i<-mod.list[[i]]
               ##Works on maxLik, lm, glm, game, and surv
               if("maxLik" %in% class(mod.i)){ ##Check if the class is maxLik
                    out.L<-str_c(out.L,
                                 paste(" & ",
                                       "\\multicolumn{1}{c}{",
                                       num2str(logLik(mod.i)),
                                       "}"))
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R, " & ")
                    out.N<- str_c(out.N,
                                  paste(" & ",
                                        "\\multicolumn{1}{c}{", 
                                        nrow(mod.i$gradientObs),
                                        "}"))
               }
               if("lm" == class(mod.i)[1]){
                    out.L<-str_c(out.L, " & ")
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R,
                                 paste(" & ",
                                       "\\multicolumn{1}{c}{", 
                                       num2str(summary(mod.i)$adj.r.squared)),
                                 "}")
                    out.N<-str_c(out.N,
                                 paste("& ", 
                                       "\\multicolumn{1}{c}{",
                                       nrow(mod.i$model),
                                       "}"))
               }
               if("plm" == class(mod.i)[1]){
                    out.L<-str_c(out.L, " & ")
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R,
                                 paste(" & ",
                                       "\\multicolumn{1}{c}{", 
                                       num2str(summary(baseFETE)$r.squared[2]),
                                 "}"))
                    out.N<-str_c(out.N,
                                 paste("& ", 
                                       "\\multicolumn{1}{c}{",
                                       nrow(mod.i$model),
                                       "}"))
               }
               if("ivreg" == class(mod.i)[1] ){
                    out.L<-str_c(out.L, " & ")
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R, " & ")
                    out.N<- str_c(out.N,
                                  paste(" & ", 
                                        "\\multicolumn{1}{c}{",
                                        nrow(mod.i$model),
                                        "}"))
               }
               
               if("glm" == class(mod.i)[1]  | "game" %in% class(mod.i) 
                  |"polr" %in% class(mod.i)){
                    out.L<-str_c(out.L,
                                 paste(" & ",
                                       "\\multicolumn{1}{c}{",
                                       num2str(as.numeric(logLik(mod.i))),
                                       "}"))
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R, " & ")
                    out.N<- str_c(out.N,
                                  paste(" & ",   
                                        "\\multicolumn{1}{c}{",
                                        nrow(mod.i$model),
                                        "}"))
               }
               if("negbin" == class(mod.i)[1]){
                    out.L<-str_c(out.L,
                                 paste(" & ",
                                       "\\multicolumn{1}{c}{",
                                       num2str(mod.i$twologlik/2),
                                       "}"))                
                    out.A <- str_c(out.A, 
                                   paste(" & ",    
                                         "\\multicolumn{1}{c}{",
                                         num2str(1/mod.i$theta),
                                         "}"))
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R, " & ")
                    out.N<- str_c(out.N,
                                  paste(" & ",
                                        "\\multicolumn{1}{c}{", 
                                        nrow(mod.i$model),
                                        "}"))
               }
               if("polywog" == class(mod.i)[1]){
                    out.L<-str_c(out.L, " & ")
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, 
                                   paste(" & ",                                       
                                         "\\multicolumn{1}{c}{", 
                                         num2str(mod.i$lambda),
                                         "}"))
                    out.Rho<-str_c(out.Rho, " & ")
                    out.R<-str_c(out.R, " & ")
                    out.N<- str_c(out.N,
                                  paste(" & ",
                                        "\\multicolumn{1}{c}{", 
                                        mod.i$nobs,
                                        "}"))
               }
               if("spautolm" == class(mod.i)[1]){
                    out.L<-str_c(out.L,
                                 paste(" & ",
                                       "\\multicolumn{1}{c}{",
                                       num2str(mod.i$LL),
                                       "}"))
                    out.A<-str_c(out.A, " & ")
                    out.Lam<-str_c(out.Lam, " & ")
                    out.Rho<-str_c(out.Rho, 
                                   paste(" & ",                                       
                                         "\\multicolumn{1}{c}{", 
                                         num2str(mod.i$lambda),
                                         "}"))
                    out.R<-str_c(out.R, " & ")
                    out.N<- str_c(out.N,
                                  paste(" & ",                                       
                                        "\\multicolumn{1}{c}{", 
                                        nrow(mod.i$X),
                                        "}"))
               }
               if("list" == class(mod.i)[1]){
                    if(is.null(mod.i$logLik)){
                         out.L<-str_c(out.L, " & ")
                         out.A<-str_c(out.A, " & ")
                         out.Lam<-str_c(out.Lam, " & ")
                         out.Rho<-str_c(out.Rho, " & ")
                         out.R<-str_c(out.R,  paste(" & ",                                       
                                                    "\\multicolumn{1}{c}{",
                                                    num2str(mod.i$ar2),
                                                    "}"))
                         out.N<- str_c(out.N,
                                       paste(" & ",                                       
                                             "\\multicolumn{1}{c}{", 
                                             mod.i$Obs,
                                             "}"))
                    }else{
                         out.L<-str_c(out.L,
                                      paste(" & ",                                       
                                            "\\multicolumn{1}{c}{",
                                            num2str(mod.i$logLik),
                                            "}"))
                         out.A<-str_c(out.A, " & ")
                         out.Lam<-str_c(out.Lam, " & ")
                         out.Rho<-str_c(out.Rho, " & ")
                         out.R<-str_c(out.R, " & ")
                         out.N<- str_c(out.N,
                                       paste(" & ",                                       
                                             "\\multicolumn{1}{c}{", 
                                             mod.i$Obs,
                                             "}"))
                    }
               }
               modNames <- str_c(modNames, "& ", model.names[i])
               modNums <- str_c(modNums, "& ", model.numbers[i])
               depNames <- str_c(depNames, " & ", dep.varnames[i])
          }
          outAll <- c(out.L, out.A, out.Lam, out.Rho, out.R, out.N)
          outI <- c(str_detect(out.L, "multi"),
                    str_detect(out.A, "multi"),
                    str_detect(out.Lam, "multi"),
                    str_detect(out.Rho, "multi"),
                    str_detect(out.R, "multi"),
                    str_detect(out.N, "multi"))
                 
          out <- paste(outAll[which(outI)], "\\\\")
          out[1] <- paste("\\hline", out[1])
          out[length(out)] <- paste(out[length(out)], "\\hline")
          out <- str_c(out, collapse="\n")
      
          stars <- ifelse(stars=="all",
                          paste("\\multicolumn{", 
                                1+length(mod.list),
                                "}{l}{\\footnotesize{\\emph{Notes:}  $^{***}p<0.01$; $^{**}p<0.05$; $^{*} p<0.10$ }} \\\\ \\multicolumn{",
                                1+length(mod.list),
                                "}{l}{\\footnotesize{", 
                                seNote,
                                "}}%",
                                sep=""),
                          ifelse(stars=="default",
                                 paste("\\multicolumn{", 
                                       1+length(mod.list),
                                       "}{l}{\\footnotesize{\\emph{Notes:} $^{*}p<0.05$}}\\\\ \\multicolumn{",
                                       1+length(mod.list),
                                       "}{l}{\\footnotesize{", 
                                       seNote,
                                       "}}%",
                                       sep=""),
                                 paste("\\multicolumn{", 
                                       1+length(mod.list),
                                       "}{l}{\\footnotesize{\\emph{Notes:} $^{***}$ 99\\% CI does not contain 0; $^{**}$ 95\\% CI does not contain 0; $^{*}$ 90\\% CI does not contain 0}}\\\\ \\multicolumn{",
                                       1+length(mod.list),
                                       "}{l}{\\footnotesize{", 
                                       seNote,
                                       "}}%",
                                       sep="")))
          
          
          header <- str_c(depNames, " \\\\ \n ", 
                          modNames, " \\\\ \n",
                          modNums, " \\\\ \n")
          
          out<-list(header=header,
                    out=str_c(out, stars))
     }
     if(is.null(out)){stop("Object class not supported")}
     return(out)
}

##handrolled Coeftest
coeftest2<-function(beta, sd){
     z<-beta/sd
     p.value<-pnorm(abs(z), lower.tail=FALSE)*2
     results<-cbind(beta, sd,z, p.value)
     colnames(results)<-c("Estimate",
                          "Std. Error",
                          "Z Value",
                          "Pr(|z|)")
     return(results)
}

##A wrapper for xtable to allow stackign multiple models and inputing 
##Bootstrapped standard errors
xtable2<-function(mod.list, 
                   se=NULL, 
                   info.list=NULL, #Only if coef=TRUE
                   coef=FALSE,
                   apsr=FALSE, 
                   bootCI=FALSE,
                   stars=c("default", "all", "CI"),
                   caption=NULL,
                   label=NULL, 
                   align=NULL, 
                   digits=2,
                   order=NULL,
                   seNote="Standard Errors in Parenthesis",
                   covariate.labels=NULL,
                   model.numbers=NULL,
                   model.names=NULL,
                   dep.varnames=NULL,         
                   k=1,
                   print.xtable.options=list()
                   ){
     
     if(length(stars)>1){
          stars <- stars[1]
     }
     
     
     if(coef){
          if(length(dep.varnames)==length(info.list)){
               dep.varnames <- paste("\\multicolumn{1}{c}{", dep.varnames, "}")
          }
          if(length(model.names)==length(info.list)){
               model.names <- paste("\\multicolumn{1}{c}{", model.names, "}")
          }
          if(length(model.numbers)==length(info.list)){
               model.numbers <- paste("\\multicolumn{1}{c}{", model.numbers, "}")
          }
               
          info <- mod.sum(info.list, apsr, stars, seNote, model.numbers, model.names, dep.varnames,k)
          if(bootCI){
            info$length <- c(0, 2*length(unique(unlist(lapply(mod.list, 
                                                              function(x){colnames(x)})))))
          }else{
            info$length <- c(0, 2*length(unique(unlist(lapply(mod.list, 
                                                              function(x){names(x)})))))
          }
          
     }else{
          if(length(dep.varnames)==length(mod.list)){
               dep.varnames <- paste("\\multicolumn{1}{c}{", dep.varnames, "}")
          }
          if(length(model.names)==length(mod.list)){
               model.names <- paste("\\multicolumn{1}{c}{", model.names, "}")
          }
          if(length(model.numbers)==length(mod.list)){
               model.numbers <- paste("\\multicolumn{1}{c}{", model.numbers, "}")
          }
          info <- mod.sum(mod.list, apsr, stars, seNote, model.numbers, model.names, dep.varnames,k)
          info$length <- c(0, 2*length(unique(unlist(lapply(mod.list, 
                                                            function(x){names(coef(x))})))))
          
     }
     
     if(apsr){
          print.xtable.default <- list(booktabs=TRUE,
                                       sanitize.text.function=function(x){x},
                                       include.rownames=FALSE,
                                       caption.placement="top",
                                       table.placement="h!",
                                       include.colnames =FALSE,
                                       add.to.row=list(pos=list(info[[3]][1],info[[3]][2]),
                                                       command=c(info[[1]], info[[2]]))
                                       )
     }else{
          print.xtable.default <-list(booktabs=TRUE,
                                      sanitize.text.function=function(x){x},
                                      include.rownames=FALSE,
                                      caption.placement="top",
                                      table.placement="h!",
                                      include.colnames =FALSE,
                                      add.to.row=list(pos=list(info[[2]]),
                                                      command=info[[1]])
          )
          
     }
     
     print.xtable.options <- modifyList(print.xtable.default, print.xtable.options)
     
     
     
     if(!coef){
          coef.list<-lapply(mod.list, coef)
     }else{
          if(bootCI){
               coef.list<-lapply(mod.list, apply, 2, mean, na.rm=TRUE)
               }else{
                coef.list<-mod.list
          }
     }
     if(is.null(order)){
         order <- unique(unlist(lapply(coef.list, names)))

     }else{
          order <- unique(c(order, unlist(lapply(coef.list, names))))
     }
     order <- rep(order, each=2)
     order[seq(2, length(order), by=2)] <- paste(order[seq(2, length(order), by=2)], 
                                                 "_SE", 
                                                 sep="")
     
     output<-data.frame()
     for(i in 1:nrow(summary(mod.list))){
          
          if(!is.null(se)){
               table<-coeftest2(coef.list[[i]], se[[i]])
          }else{
               if(bootCI){
                    table <- cbind(coef.list[[i]], 
                                   t(apply(mod.list[[i]],
                                           2, 
                                           quantile,
                                           c(0.005, 0.995),
                                           na.rm=TRUE)),
                                   t(apply(mod.list[[i]],
                                           2, 
                                           quantile,
                                           c(0.025, 0.975),
                                           na.rm=TRUE)),
                                   t(apply(mod.list[[i]],
                                           2, 
                                           quantile,
                                           c(0.05, 0.95),
                                           na.rm=TRUE)))
               }else{
                    table<-coeftest2(coef.list[[i]], 
                                     sqrt(diag(vcov(mod.list[[i]]))))
               }
          }
          
          if(apsr){
               a.table<-data.frame()
               for(j in 1:nrow(table)){
                    if(!bootCI){
                         temp<-table[j, 1:2]
                         temp<-t(t(temp))
                         temp<-num2str(temp)
                         temp<-str_trim(temp)
                         temp[2,]<-paste("(",
                                    temp[2,],
                                    ")", sep="")
                    }else{
                         temp <- table[j, c(1,4:5)]
                         temp <- num2str(temp)
                         temp[2] <- paste("\\myalign{,}{(",
                                          str_c(temp[2:3], 
                                                collapse=", "),
                                          ")}",
                                          sep="")
                         temp <- temp[-3]
                         temp <- t(t(temp))
                    }
                    if(stars=="default"){
                         if(table[j, 4]<0.05){
                              temp[1,]<- paste(temp[1,] ,"$^*$", sep="")
                         }
                    }else{
                         if(stars=="all"){
                              if(table[j, 4]<0.1 & table[j, 4]>0.05){
                                   temp[1,]<- paste(temp[1,] 
                                                    ,"$^{*}$", sep="")
                              }
                              if(table[j, 4]<0.05 & table[j, 4]>0.01){
                                   temp[1,]<- paste(temp[1,],
                                                    "$^{**}$", sep="")
                              }
                              if(table[j, 4]<0.01){
                                   temp[1,]<- paste(temp[1,],
                                                    "$^{***}$", sep="")
                              }
                         }else{
                              if(sign(table[j,2])==sign(table[j,3])){
                                   temp[1,] <- paste(temp[1,],
                                                     "^{***}",
                                                     sep="")
                              }else{
                                   if(sign(table[j,4])==sign(table[j,5])){
                                        temp[1,] <- paste(temp[1,],
                                                          "^{**}", 
                                                          sep="")
                                   }else{
                                        if(sign(table[j,6])==sign(table[j,7])){
                                             temp[1,] <- paste(temp[1,],
                                                               "^{*}", 
                                                               sep="")
                                        }
                                   }
                              }
                         }
                    }
                    row.names(temp)<-NULL
                    var.name<-row.names(table)[j]
                    colnames(temp)<-paste("Model.", i, sep="")
                    row.names(temp)<-c(var.name, paste(var.name, "_SE", sep=""))
                    a.table<-rbind(a.table, temp)  
               }
               if(i==1){
                    a.output<-a.table
                    ` `<-row.names(a.output)
                    row.names(a.output)<-NULL
                    a.output<-cbind(` `, a.output)
                    if(is.null(order)){
                         a.output$` ` <- factor(a.output$` `, 
                                               levels=as.character(a.output$` `))
                    }else{
                         a.output$` ` <- factor(a.output$` `, 
                                               levels=order)
                    }
               }else{
                    ` `<-row.names(a.table)
                    row.names(a.table)<-NULL
                    a.table<-cbind(` `, a.table) 
                    if(is.null(order)){
                         a.table$` ` <- factor(a.table$` `, 
                                               levels=as.character(a.table$` `))
                    }else{
                         a.table$` ` <- factor(a.table$` `, 
                                               levels=order)
                    }
                    
                    a.output<-merge(a.output, 
                                    a.table,
                                    by= " ",
                                    all=TRUE,
                                    sort=TRUE)        
                    a.output <- a.output[order(a.output[,1]),]
               }
          }else{
               ` `<-row.names(table)
               row.names(table)<-NULL
               table<-cbind.data.frame(` `, table)
               output<-rbind(output, table)
               int <- output[1,]
               output <- output[2:nrow(output),]
               output <- rbind(output, int)
          }
     }
     if(apsr){
          a.output$` ` <- as.character(a.output$` `)
          a.output$` `[str_detect(a.output$` `, "_SE")]<-""
          output<-a.output
          const <- str_detect(output[,1], ignore.case("intercept")) | 
                     str_detect(output[,1], ignore.case("constant")) |
                     str_detect(output[,1], ignore.case("const"))
            
          const[is.na(const)] <- FALSE
          idx <- which(const)
          int <- output[c(idx, idx+1),]
          output <- output[-c(idx, idx+1),]
          output <- rbind(output, int)
          colnames(output) <- rep("", ncol(output))    
     }
     if(class(covariate.labels)=="character"){
          output[output[,1]!="",1] <- covariate.labels
     }
     align<-ifelse(is.null(align),
                   ifelse(apsr==FALSE, 
                          str_c(rep("r", ncol(output)+1), collapse=""),
                          ifelse(!bootCI,
                                 str_c(c("rr", rep("S", ncol(output)-1)), collapse=""),
                                 str_c(c("rr", rep("d", ncol(output)-1)), collapse=""))),
                   align)
     #      return(output)
     suppressWarnings(
     x <- xtable(output,
                 caption=caption,
                 label=label, 
                 align=align, 
                 digits=digits))
  print.xtable.options <- modifyList(print.xtable.options, list(x=x))
  do.call(print, print.xtable.options)
     
     if(bootCI){
          cat("Reminder: Add the following to your TeX preamble:\n")
          cat("\\usepackage{dcolumn}\n
\\newcolumntype{d}{D{.}{.}{-1}}\n
\\newcolumntype{,}{D{,}{,\\,}{-1}}\n
\\newcommand*{\\myalign}[2]{\\multicolumn{1}{#1}{#2}}\n")
     }
     if(str_detect(align, "S")){
          cat("Reminder: Add the following to your TeX preamble:\n")
          cat("\\usepackage{siunitx}
\\sisetup{\n
     input-symbols=(),\n
	table-align-text-post = false,\n
	group-digits=false,\n
} ")
     }
}
###Note this function requires 


#\usepackage{siunitx}
#Be in the preamble 



av.data<-function(x, x1=NULL, length.out=25){ ##New function to set columns in
     ##a matrix or data.frame at their mean (for continuous variables)
     ##or median (dummies)
     ##ARGUMENTS:
     ##x=Data
     ##x1= variable that should vary (if any)
     ##length.out=nrow of output
     
     
     if(class(x)=="data.frame"){
          x<-as.matrix(x) 
     }##convert data.frames into matrices (I just like them better)
     
     newdata<-x ##want dim(newdata)==dim(x)
     for(i in 1:ncol(x)){
          if(all(x[,i] %in% c(0,1))){ ##find dummies
               newdata[,i]<-median(x[,i]) ##set to median
          }else{ ##not a dummy?
               newdata[,i]<-mean(x[,i]) ##Set to mean
          }     
     }
     
     newdata<-newdata[1:length.out,] ##Trim to desired length
     
     if(!is.null(x1)){ ##if user inputs an x1, then
          X1<-seq(min(x[,x1]),  ##Create a sequence from min
                  max(x[,x1]),  ##to max
                  length.out=length.out) ##with desired length
          newdata[,x1]<-X1 ##And put it in the right column
     }
     
     
     
     return(newdata)
}

#######Likelihood Functions######
NBreg<-function(b, X, Y, single=FALSE){
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     
     alpha<-exp(b[length(b)])
     beta<-b[-length(b)]
     e.XB<-exp(X %*% beta)
     LL<-Y*log((alpha*e.XB)/(1+alpha*e.XB))-(1/alpha)*log(1+alpha*e.XB) + 
          lgamma(Y + 1/alpha)-lgamma(Y+1)-lgamma(1/alpha)
     if(single){
          return(sum(LL))
     }else{
          return(LL)
     }
}

Pois <- function(b, X, Y){
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     LL <- Y * (X%*%b) - exp(X %*% b)
     return(LL)
}

Probit<-function(b, X, Y){
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     Xb<-X %*% b
     ll<-ifelse(Y==1,
                pnorm(Xb, log.p=TRUE),
                pnorm(-Xb, log.p=TRUE))
     return(ll)
} ##Probit likelihood function

Llogit=function(b,X,Y){
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     xb=X%*%b
     logD=ifelse(-xb>708,-xb,log(1+exp(-xb)))
     llik=(1-Y)*(-xb)-logD
     return(llik)
}
##FE neg bin
FEnbreg <- function(b, X, Y, FE){
     if(!all(X[,1]==1)){X <- cbind(1, X)}
     llik <- matrix(ncol=1, nrow=0)
     XB <- X %*% b 
     for(i in unique(FE)){
          llik <-  rbind(llik,
                             lgamma(sum(exp(XB[FE==i])))  + 
                                  lgamma(sum(Y[FE==i])+1)- 
                                 # sum(lgamma(Y[FE==i]+1))- 
                                  lgamma(sum(Y[FE==i])+ sum(exp(XB[FE==i])))+
                                  sum(lgamma(exp(XB[FE==i])+Y[FE==i]) -
                                  lgamma(exp(XB[FE==i]))- lgamma(Y[FE==i]+1)))
     }
     return(llik)
}


grFEnbreg <- function(b, X, Y, FE){
     if(!all(X[,1]==1)){X <- cbind(1, X)}
     gr <- matrix(ncol=ncol(X), nrow=0)
     XB <- X %*% b 
     for(i in unique(FE)){
          gr <-  rbind(gr,
                         colSums(X[FE==i,,drop=FALSE]* exp(XB[FE==i])) * digamma(sum(exp(XB[FE==i])))  - 
                         colSums(X[FE==i,,drop=FALSE]* exp(XB[FE==i])) * digamma(sum(Y[FE==i])+ sum(exp(XB[FE==i])))+
                         colSums(X[FE==i,,drop=FALSE]* exp(XB[FE==i]) *  digamma(exp(XB[FE==i])+Y[FE==i])  -
                         X[FE==i,,drop=FALSE]* exp(XB[FE==i]) *  digamma(exp(XB[FE==i])) ) )
     }
     return(gr)
}

diapprox <- function(x){
     return((lgamma(x + 0.0001)- lgamma(x-0.0001))/0.0002)
}



NBgrad <- function(b, X, Y){ ##DO NOT USE
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     
     eta <- b[length(b)]
     alpha<-exp(eta)
     a2 <- exp(2*eta)
     beta<-b[-length(b)]
     XB <- as.numeric(X %*% beta)
     e.XB<-exp(XB)
     e.2XB<-exp(2*XB)
     eXBA <- exp(XB+eta)
     
    dB <- X* (Y-e.XB)/(eXBA+1)
    dA <-( exp(-eta)/(eXBA+1)) * (Y * alpha + 
                                      eXBA * log(eXBA+1) - 
                                      eXBA * digamma(Y + 1/alpha) +
                                      eXBA * digamma(1/alpha) -
                                      eXBA+
                                      log(eXBA + 1)-
                                      digamma(Y + 1/alpha) + 
                                      digamma(1/alpha))
    return(cbind(dB,dA))
                                      
}

##Weighted relogit
weighted <- function(b, X, Y, w){
     if(!all(X[,1]==1)){X <- cbind(1, X)}
     return(-w * log(1 + exp((1-2*Y) * (X %*%b))))
}



########Comparative Model Testing --Make sure these work before using######
Clarke<-function(mod1, mod2, lik1, lik2){
     if("maxLik" %in% class(mod1) & 
             "maxLik" %in% class(mod2)){
          n<-nrow(mod1$gradientObs)     
          test<-(lik1-length(mod1$estimate)*log(n)/(n*2) )>
               (lik2-length(mod2$estimate)*log(n)/(n*2))
          p<-pbinom(sum(test),size=n, p=0.5, lower.tail=F)
          results<-matrix(c(sum(test), p), nrow=1)
          colnames(results)<-c("Test Statistic",
                               "$p$-value")
          row.names(results)<-"Clarke test results"
          return(results)
     }else{stop("only maxLik supported")}
}

Vuong<-function(mod1, mod2, lik1, lik2){
     if("maxLik" %in% class(mod1) & 
             "maxLik" %in% class(mod2)){
          LR<-logLik(mod1)-logLik(mod2)
          dif<-lik1-lik2
          n<-nrow(mod1$gradientObs)
          SC<-(length(mod1$estimate)- 
                    length(mod2$estimate))*log(n)/2
          w<-sqrt(mean(dif^2)-mean(dif)^2)*sqrt(n)
          test<-(LR-SC)/w
          test<-(LR)/w
          p<-2*pnorm(abs(test), lower.tail=FALSE)
          results<-matrix(c(LR, test, p), nrow=1)
          colnames(results)<-c("Likelihood Ratio",
                               "Test Statistic",
                               "$p$-value")
          row.names(results)<-"Vuong test results"
          return(results)
     }else{stop("only maxLik supported")}
}



###########Moran's I test################
MoranI <- function(u, W){
     N <- nrow(W)
     S1 = .5 * sum((W + t(W))^2)
     S2 = sum((colSums(W) + rowSums(t(W)))^2)
     S3 = ((1/49) * sum(u -mean(u))^4)/((1/49) * sum((u - mean(u))^2))^2
     S4 = (N^2 - 3*N +3)*S1 - N*S2 + 3*(sum(W)^2)
     S5 = S1 - 2*N*S1 + 6*(sum(W)^2)
     
     V =  (N*S4 - S3*S5)/( (N-1 )*(N-2) *(N-3) * (sum(W)^2))
     I =(N / sum(W)) * sum(W* (u-mean(u)) %o% (u-mean(u))) / sum((u-mean(u))^2)
     Ei <- -1/(N-1)
     
     test <- abs(I- Ei)/sqrt(V)
     p <- (1-pnorm(test, lower.tail=TRUE))*2
     
     results <- list(Observed = I, Expected=Ei, s.e.= sqrt(V), Z= test, p =p )
     return(results)
}


#####NegBin Clarke Test#######
clarke.nb <- function(x,z){
     x <- nb4
     z <- nb2
     
     y <- x$y
     N <- length(y)
     a1 <- 1/x$theta 
     a2 <- 1/z$theta
     fit.mod.1 <- fitted(x)
     fit.mod.2 <- fitted(z)
     
     xb.1 <- predict.glm(x)
     xb.2 <- predict.glm(z)
     
     
     
     nbll <- function(xb, Y, alpha){
          e.XB<-exp(xb)
          LL<-Y*log((alpha*e.XB)/(1+alpha*e.XB))-(1/alpha)*log(1+alpha*e.XB) + 
               lgamma(Y + 1/alpha)-lgamma(Y+1)-lgamma(1/alpha)
          return(LL)
     }
     
     loglik.1 <- nbll(xb.1, y, a1)
     loglik.2 <- nbll(xb.2, y, a2)
     
     pdim <- x$rank
     qdim <- z$rank
     adjll.1 <- loglik.1-(((pdim)/(2*N))*log(N))
     adjll.2 <- loglik.2-(((qdim)/(2*N))*log(N))
     adjdiff <- ((adjll.1-adjll.2)>0) 
     
     Clarke.test <- binom.test(sum(adjdiff),N,p=.5)
     Clarke <- Clarke.test$statistic
     Clarke.p <- Clarke.test$p.value

     ClarkeAns <- list(Statistic = Clarke,
                       p.value = Clarke.p)
     return(ClarkeAns)
}




# #### Part of modinfo that I don't think I need anymore
# if(all(sapply(lapply(lapply(mod.list, class), `%in%`, c("lm", "plm", "panelmodel")), all))){
#      out<-str_c(out.R, "\\\\\n", 
#                 out.N,   "\\\\ \\hline  \n")
# }else{
#      ifelse(any("negbin" %in% unlist(lapply(mod.list, class))),
#             ifelse(any("spautolm" %in% unlist(lapply(mod.list, class))),
#                    ifelse(any("polywog" %in% unlist(lapply(mod.list, class))),
#                           ifelse(any(c("plm", "lm")%in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.R, "\\\\\n", 
#                                             out.Lam, "\\\\\n", 
#                                             out.Rho, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.Lam, "\\\\\n", 
#                                             out.Rho, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n")),
#                           ifelse(any(c("plm", "lm") %in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.Rho, "\\\\\n",
#                                             out.R, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.Rho, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n"))),
#                    ifelse(any("polywog" %in% unlist(lapply(mod.list, class))),
#                           ifelse(any(c("plm", "lm") %in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.R, "\\\\\n", 
#                                             out.Lam, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.Lam, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n")),
#                           ifelse(any(c("lm", "plm") %in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.R, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.A, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n")))),
#             ifelse(any("spautolm" %in% unlist(lapply(mod.list, class))),
#                    ifelse(any("polywog" %in% unlist(lapply(mod.list, class))),
#                           ifelse(any(c("plm", "lm") %in%unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.R, "\\\\\n", 
#                                             out.Lam, "\\\\\n", 
#                                             out.Rho, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.Lam, "\\\\\n", 
#                                             out.Rho, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n")),
#                           ifelse(any(c("plm", "lm") %in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.Rho, "\\\\\n",
#                                             out.R, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.Rho, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n"))),
#                    ifelse(any("polywog" %in% unlist(lapply(mod.list, class))),
#                           ifelse(any(c("plm", "lm") %in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.R, "\\\\\n", 
#                                             out.Lam, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.Lam, "\\\\\n",
#                                             out.N,   "\\\\ \\hline  \n")),
#                           ifelse(any(c("plm", "lm") %in% unlist(lapply(mod.list, class))),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.R, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n"),
#                                  out<-str_c(out.L, "\\\\\n", 
#                                             out.N,   "\\\\ \\hline  \n")))))
# } 
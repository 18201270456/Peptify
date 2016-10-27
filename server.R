library(shiny)
attach("myWorkspace.RData")
library(KernSmooth)

shinyServer(function(input, output, session) {
  
  inputData <- reactive({  
      input$inputId
  })
  
  inputBoxTrain <- reactive({  
    input$train
  })
  
  inputAA <- reactive({  
    input$aa
  })
  
  inputBG <- reactive({  
    input$bg
  })
  
  output$plot <- renderPlot({
    #input1 = gsub(" ","",inputData()) 
    if(!is.null(inputData())) peptides = strsplit(inputData(), split = " ") else peptides = ""
    peptideClean = as.vector(toupper(peptides))
    # cat(inputData()[1])
    # cat(inputData()[2])
    # cat(inputData())
    # cat(inputBoxTrain()[1])
    try(
     peptify(peptideClean, showtrainingset=inputBoxTrain(),showaa=inputAA(), showbackground=inputBG(), showplot='true')
    ,silent = TRUE
    )
  })
  
})

peptify <- function(peptide, showplot='true', showtrainingset='true', showaa='false', showbackground='true') {
  if(peptide == ""){
   whitePlaceholder = TRUE
  } else {
    whitePlaceholder = FALSE
    
    limits=c(-1,1)
    PREDICTED_SPLITTED = c()
    PREDICTED_SPLITTED = strsplit(peptide,'')
    PP1_PREDICTED=c()
    PP2_PREDICTED=c()
    CPP_PREDICTED=c()
    
    for (i in 1:length(PREDICTED_SPLITTED)) {
      PP1_PREDICTED[i]=mean(DESC[PREDICTED_SPLITTED[[i]][],'PP1'])
      PP2_PREDICTED[i]=mean(DESC[PREDICTED_SPLITTED[[i]][],'PP2'])
    }
    
  }
  
  if (showplot) {
      
      plot(PP2_CPP, PP1_CPP,type='n', xlim=c(-1,1), ylim=c(-1,1), xlab="(Hydrophobicity)", ylab="(Polarity)")
      
      if (showbackground) {
        xx=c(-1,-0.764, 0.41,-1)
        yy=c(-1,-1, 1,1)
        polygon(xx,yy,col="#F7DDDC", border=NA)
        
        xx=c(-0.764,1, 1,0.41)
        yy=c(-1,-1, 1,1)
        polygon(xx,yy,col="#F5EAB3", border=NA)
        
        xx=c(0,1, 1,0)
        yy=c(-0.3,-0.3, 1,1)
        polygon(xx,yy,col="#D7F2BD", border=NA)
        
        xx=c(-0.47,1, 1,0.71)
        yy=c(-1,-1, 1,1)
        polygon(xx,yy,col="#D7F2BD", border=NA)
      }
      
      if (showtrainingset) {
        points(PP2_DECOY[1:5000], PP1_DECOY[1:5000], col = densCols(PP2_DECOY[1:5000], PP1_DECOY[1:5000], nbin = 100, colramp=colorRampPalette(c("#e6b5b5","#610a0a"))), pch=20)
        points(PP2_CPP, PP1_CPP, col='#247A0F', pch=20)
      }
      
      if (showaa) {
        text(DESC$PP2,DESC$PP1,DESC$AA1)
      }
  }

  if(whitePlaceholder == FALSE){  
    
    for (i in 1:length(PREDICTED_SPLITTED)) {
      if ((PP2_PREDICTED[i]*1.7 + 0.3)>PP1_PREDICTED[i]) {
        CPP_PREDICTED[i]=0
      } else {
        CPP_PREDICTED[i]=-1
      }
      if ((PP2_PREDICTED[i]*1.7 - 0.2)>PP1_PREDICTED[i] | PP2_PREDICTED[i]>0) {
        CPP_PREDICTED[i]=1
      } 
    }
    
    points(PP2_PREDICTED, PP1_PREDICTED, col='red', pch=18, cex = 1.8)
    
    }
}

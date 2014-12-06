library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  output$plotout <- renderPlot({
    td <- make.treedata(tree, filter_(dat, paste("TEST == '", input$test,"'", sep=""))) 
    td <- filter(td, !outliers(ISIS_MEAN, input$iqrfactor))
    if(input$logdat){
      logfn <- log
      yaxis <- paste("Log", input$test, sep=" ")
    } else {
      logfn <- function(x) x
      yaxis <- input$test
    }
    y <- logfn(td$dat$ISIS_MEAN)
    phyhalf <- log(2)/phylolm(y~1,phy=td$phy , model="OUfixedRoot")$optpar
    if(input$plotfn=="Phenogram"){
      phenogram(td$phy, logfn(setNames(td$dat$ISIS_MEAN, attributes(td)$tip.label)), fsize=0.5, ylab=yaxis,main=paste("phylogenetic half-life:", round(phyhalf,4)))
    }
    if(input$plotfn=="Ancestral States"){
      td1 <- select(td, ISIS_MEAN)
      td1$dat[[1]] <- logfn(td1$dat[[1]])
      res <- aceArbor(td1, charType="continuous", aceType="marginal")
      TH <- max(branching.times(td$phy))
      plot(res, label.offset=0.05*TH, cex.asr=1, cex=0.5, main=paste("phylogenetic half-life:", round(phyhalf,4)))
    }
      })
})
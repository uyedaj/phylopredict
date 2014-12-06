library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Normal animal reference values"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("test",
                  "TEST:",
                  as.character(unique(dat$TEST))),
      selectInput("logdat",
                  "LOG:",
                  c(FALSE, TRUE)),
      selectInput("plotfn",
                  "Plot Type:",
                  c("Phenogram", "Ancestral States")),
      sliderInput("iqrfactor",
                  "Exclude outliers X times outside IQR:",
                  min=1.5,
                  max=10, 
                  value=10,
                  step=1)
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plotout",height="600px")
    )
  )
))
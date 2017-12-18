
#libraries
library(shiny)
library(pwr)
library(ggplot2)
library(knitr)
library(DT)

# general functions for use
cohen_d <- function(d1,d2) {  
  m1 <- mean(d1, na.rm=TRUE)
  m2 <- mean(d2, na.rm=TRUE)
  s1 <- sd(d1, na.rm=TRUE)
  s2 <- sd(d2, na.rm=TRUE)
  spo <- sqrt((s1**2 + s2**2)/2)
  d <- (m1 - m2)/spo
  effsi <- d / sqrt((d**2)+4)
  ret <- list("d" = d, "effectsi" = effsi)
  return(ret)  

  } 

# right now this only works with the typical alpha levels
# I don't think we need to allow more precision than that but conceivable 
# people might want to look at alpha levels in between those typical thresholds.

# look at how 


dat_MDE <- function(mean.input, sd.input, differs){
  #initialise empty vec
  p <- matrix(NA, nrow = 20, ncol=3)
  
  set.seed(20171101)
  for(i in seq(1:length(differs)) ) {
    samp1 <- rnorm(n=1000, mean = mean.input, sd=sd.input)
    #this is a better version if you can understand it:
    samp2 <- samp1 + rnorm(length(samp1), differs[i], differs[i]/10) #add some noise
    inp <- cohen_d(samp1, samp2)
    
    p[i,1] <- pwr.2p.test(h=inp$effectsi , sig.level=0.01, power=0.8, n=NULL)$n
    p[i,2] <- pwr.2p.test(h=inp$effectsi , sig.level=0.05, power=0.8, n=NULL)$n
    p[i,3] <- pwr.2p.test(h=inp$effectsi , sig.level=0.1, power=0.8, n=NULL)$n
    
    
  }
  
  p <- p[!is.na(p[,1]),]
  return(p)

}


# 
# 
# mde_tab = data.frame("MDE"=differs, 
#                      "Sample size 0.05"=p[,1], 
#                      "Sample size 0.1" = p[,2]
# )
# return(mde_tab)

# to get measurable effect sizes
posthoc_mde = function(n_length){
  require(pwr)
  
  #print(paste("mean diff (d2-d1)", mean(d2, na.rm=TRUE)-mean(d1,na.rm=TRUE)))
  # what is happening here?
  # s1= sd(d1,na.rm=TRUE)
  # s2=sd(d2,na.rm=TRUE)
  # spo = sqrt((s1**2 + s2**2)/2)
  
  mded = pwr.t.test(n = n_length, d=NULL, sig.level = 0.05, power=0.8, type="two.sample", alternative = "two.sided")$d
  
  #mde_out = mded * spo
  return(mded)
}


# Define UI for application that draws a histogram
ui <- navbarPage("Practical Power Calculations",

   # Application title
   tabPanel("Decision Threshold",
        fluidRow(
          column(12,
            wellPanel(
              helpText("Use case: Here is where you do calulations to determine the sample required to detect a 
                       differnece of interest between treatment arms. You can do this as an aboslute
                       (magnitude) change or as a percentage change. What difference do we need to see to scale this trial? The app will automatically show you
                       sample sizes for differences greater and less than the desired differnece to create
                       a menu of options for decision making."),
              checkboxInput(inputId = "percentage", label = "Percentage change?", value = FALSE)
            )),
          column(4,
                 wellPanel(
              conditionalPanel(
                condition = "input.percentage == false",
                numericInput("sc_mean",
                             "Outcome average",
                             value = 100),
                numericInput("sc_stdev",
                             "Outcome standard deviation",
                             value = 50),
                numericInput("sc_diff_m",
                             "Magnitude change over average",
                             value = 10),
                # p("What range of changes do you want to see?"),
                # sliderInput("sc_diff_spread",
                #              "Range of permutations",
                #              min=0, max=10, value=10)
                
                # figure out how to use this to reshape the graph and table
                # I'll also need to require at least one of these to be selected
                # and print a default message if one isn't.
                checkboxGroupInput("alphaLevel",
                                   "For which alpha levels?",
                                   choices = c("p < 0.01",
                                               "p < 0.05",
                                               "p < 0.1"),
                                   selected = c("p < 0.01",
                                                "p < 0.05",
                                                "p < 0.1"))
              ), #conditional panel1
              conditionalPanel(
                condition = "input.percentage == true",
                numericInput("sc_mean",
                             "Outcome average",
                             value = 100),
                numericInput("sc_stdev",
                             "Outcome standard deviation",
                             value = 50),
                numericInput("sc_diff_p",
                             "Percentage (%) change over average",
                             value = 10),
                checkboxGroupInput("alphaLevel",
                                   "For which alpha levels?",
                                   choices = c("p < 0.01",
                                               "p < 0.05",
                                               "p < 0.1"),
                                   selected = c("p < 0.01",
                                                "p < 0.05",
                                                "p < 0.1"))
              ),
              downloadButton("downloadData", "Download Result")
            ) #wellpanel layout
          ), # column close
          column(8,
               wellPanel(
                 plotOutput("sc_plot"),
                 br(),
                 dataTableOutput("sc_table")
                        )
                )
      ),
      column(12,
             wellPanel(
               htmlOutput("nrequired")
             )
      )# fluidRow
  ), #tabPanel
   tabPanel("Minimum Detectable Effect",
            fluidRow(
              column(12,
                     wellPanel(
                       helpText("Use case: Power calculations for when you have a set sample size and 
                                  you're trying to understand the type of effect you're 
                                likely to detect with your trial. This reports changes
                                in terms of standardized effect sizes. This app will automatically show
                                results for larger and smaller sample sizes to facilitate
                                decision making.")
                       )),
              column(4,
                     wellPanel(
                         numericInput("mde_n",
                                      "Available sample size",
                                      value = 100),
                         numericInput("mde_sd",
                                      "Standard deviation of outcome",
                                      value = NULL)
                     ) #wellpanel layout
              ), # column close
              column(8,
                     wellPanel(
                       plotOutput("mde_plot"),
                       dataTableOutput("mde_table")
                     )
                  )
              )
        )# tab panel close
    # new tabpanel of the PPV work?
          
)



# Define server logic required to draw a histogram

server <- function(input, output) {
  
  # for the slider to update
  # observe({
  #   val <- input$sc_diff_m
  #   # Control the value, min, max, and step.
  #   # Step size is 2 when input value is even; 1 when value is odd.
  #   updateSliderInput(session, "sc_diff_spread", value = val,
  #                     min = floor(val/2), max = val+4, step = (val+1)%%2 + 1)
  # })
  
  # the same
  changes = sort(c(seq(-10,10, by=2), 1))
  
  # need to keep this simple 
  differs <- reactive({
    if(input$percentage==TRUE){
      req(input$sc_diff_p)
      # add in permutation on the percent changes +- 10
      
      toTest = ((input$sc_diff_p + changes)/100) * input$sc_mean
      
      # then convert to percentage
      # toTest = toTest / 100
      toTest = toTest[toTest>0]
      # toTest = toTest * input$sc_mean
      # toTest
    } else {
      req(input$sc_diff_m)
      changes = seq(-10,10, by=2) # plus or minus 10 by 2
      toTest = input$sc_diff_m + changes
      toTest = toTest[toTest>0]
  }  
})

  
  # magnitude
  sampTab <- reactive({
      req(input$sc_mean, input$sc_stdev)
      dat = round(dat_MDE(input$sc_mean, input$sc_stdev, differs()),0)
    # changes to evaluate in the data
  })
  
  
  output$sc_plot <- renderPlot({
    
    dat = sampTab()
    
      ggplot() +
          geom_point(aes(x=dat[,1], y= differs()), size=3, color="green", shape=1) +
          geom_line(aes(x=dat[,1], y=differs()), size=1.2, color="green") + 
          geom_point(aes(x=dat[,2], y= differs()), size=3, color="blue", shape=1) +
          geom_line(aes(x=dat[,2], y=differs()), size=1.2, color="blue") +
          geom_point(aes(x=dat[,3], y= differs()), size=3, color="red", shape=1) +
          geom_line(aes(x=dat[,3], y=differs()), size=1.2, color="red") +
          xlab("Sample size")+ ylab("Decision Threshold") +
          ggtitle("Decision Threshold vs. Sample Size") +
          theme(plot.title = element_text(hjust = 0.5,
                                        size = 20),
              plot.subtitle = element_text(hjust=0.5))
  })
  
  output$sc_table <- renderDataTable({
    
    dat = format(sampTab(), big.mark=",")
    # changes to evaluate in the data
    dat = as.data.frame(cbind(differs(), dat))
    names(dat) = c("Differences", "Alpha is 0.01", "Alpha is 0.05", "Alpha is 0.1")
    
    datatable(dat, rownames = FALSE)
    # selection = list(mode="single",
    # seleted = which(changes==1),
    # target = "row")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sampleSize ", Sys.Date(), ".csv", sep = "")
    },
    
      content = function(file) {
        write.csv(sampTab(), file, row.names = FALSE)
      }
  )
  
  output$nrequired <- renderUI({
    if(input$percentage==FALSE){
      dat = sampTab()
      target = input$sc_diff_m
      samp_out = round(dat[3,1],0) #because 1 is the 3 element in the changes vector >> update this to be more flexible.
      
      str1 = paste("To achieve 80% power for a decision threshold of", target, "you'll need ", samp_out, "farmers per treatment arm", sep = " ")
      str2 = "This is still a work in progress - I want this to summarize for all alpha levels - Maybe let people set or select alpha levels of interest? This is currently 
      only displaying for alpha = 0.01"
      HTML(paste(str1, str2, sep = '<br/>'))
    } else {
      dat = sampTab()
      target = input$sc_diff_p
      samp_out = round(dat[3,1],0) #because 1 is the 3 element in the changes vector >> update this to be more flexible.
      
      str1 = paste("To achieve 80% power for a decision threshold of", target, "% you'll need ", format(samp_out, big.mark=","), "farmers per treatment arm", sep = " ")
      str2 = "This is still a work in progress - I want this to summarize for all alpha levels - Maybe let people set or select alpha levels of interest? This is currently 
      only displaying for alpha = 0.01"
      HTML(paste(str1, str2, sep = '<br/>'))
    }
  })
  
  # minimum detectable effect section of server
  # add in look at getting mean if standard deviation is provided
  
  mde_changes = seq(0.1, 2, by=0.1)
  
  nOps = reactive({
    req(input$mde_n)
    nOps = input$mde_n * mde_changes
  })
  
  # minimum detectable effect
  output$mde_plot <- renderPlot({
    # only changes greater than 0
    
    dat = data.frame(d = unlist(lapply(nOps(), posthoc_mde)), n = nOps())
    
    ggplot(dat, aes(x = n, y = d)) + 
      geom_line(size = 1.2) +
      geom_point(data=dat[which(mde_changes==1),], aes(x = n, y = d), size=3) + 
      labs(x = "Sample size", y = "Minimum detectable effect (d)", 
           title = "Miniumum detectable effect for given sample size") + 
      theme(plot.title = element_text(hjust = 0.5, size = 20))
      
  })
  
  output$mde_table <- renderDataTable({
    if(is.null(input$mde_sd)){
    dat = data.frame(n = nOps(), d = round(unlist(lapply(nOps(), posthoc_mde)),3))
    names(dat) = c("Sample Size", "Effect Size")
    datatable(dat, rownames = FALSE, selection = list(mode = "single", 
                                                      selected=which(mde_changes==1),
                                                      target="row"))
    } else {
      dat = data.frame(n = nOps(), 
                       d = round(unlist(lapply(nOps(), posthoc_mde)),3))
      
      dat$mu = dat$d * input$mde_sd
      names(dat) = c("Sample Size", "Effect Size", "Average difference")
      datatable(dat, rownames = FALSE, selection = list(mode = "single", 
                                                        selected=which(mde_changes==1),
                                                        target="row"))
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


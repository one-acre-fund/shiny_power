
#libraries
library(shiny)
library(pwr)
library(ggplot2)
library(knitr)
library(DT)
library(clusterPower)
library(reshape2)

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

# add in some notes to test making commits through githubs

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
  return(as.data.frame(p))

}

# test <- crtpwr.2mean(alpha = 0.01, power = 0.8,
#              n = 10, cv = 0, d = 10, varw = 3,
#              icc = 0.01, method = "taylor")


# clustered MDE for decision threshold
# we're most likely solving for m, the number of clusters per arm as the size of the cluster,
# the site, group, etc, is likely set.
dat_MDE_clus <- function(mean.input, sd.input, clustN.input,
                         icc.input, differs){
  #initialise empty vec
  p <- matrix(NA, nrow = 20, ncol=3)
  
  set.seed(20171101)
  for(i in seq(1:length(differs)) ) {
    # samp1 <- rnorm(n=1000, mean = mean.input, sd=sd.input)
    # #this is a better version if you can understand it:
    # samp2 <- samp1 + rnorm(length(samp1), differs[i], differs[i]/10) #add some noise
    # inp <- cohen_d(samp1, samp2)
    
    p[i,1] <- crtpwr.2mean(alpha = 0.01, power = 0.8,
                           n = clustN.input, cv = 0, d = differs[i], varw = sd.input,
                           icc = icc.input, method = "taylor")[[1]]
      
    p[i,2] <- crtpwr.2mean(alpha = 0.05, power = 0.8,
                           n = clustN.input, cv = 0, d = differs[i], varw = sd.input,
                           icc = icc.input, method = "taylor")[[1]]
    
    p[i,3] <- crtpwr.2mean(alpha = 0.1, power = 0.8,
                           n = clustN.input, cv = 0, d = differs[i], varw = sd.input,
                           icc = icc.input, method = "taylor")[[1]]
    
  }
  
  p <- p[!is.na(p[,1]),]
  return(as.data.frame(p))
  
}

# to get measurable effect sizes
posthoc_mde = function(n_length){
  require(pwr)
  
  #print(paste("mean diff (d2-d1)", mean(d2, na.rm=TRUE)-mean(d1,na.rm=TRUE)))
  # what is happening here?
  # s1= sd(d1,na.rm=TRUE)
  # s2=sd(d2,na.rm=TRUE)
  # spo = sqrt((s1**2 + s2**2)/2)
  spo = 2 # this is crude but need to have some sort of inflation factor until I figure out why we need this.
  
  a1 = pwr.t.test(n = n_length, d=NULL, sig.level = 0.01, power=0.8, type="two.sample", alternative = "two.sided")$d
  a5 = pwr.t.test(n = n_length, d=NULL, sig.level = 0.05, power=0.8, type="two.sample", alternative = "two.sided")$d
  a10 = pwr.t.test(n = n_length, d=NULL, sig.level = 0.1, power=0.8, type="two.sample", alternative = "two.sided")$d
  
  mded = as.data.frame(cbind(a1, a5, a10))
  
  mded = mded * spo
  return(mded)
}

# start UI function
ui <- navbarPage("Practical Power Calculations",

   # Application title
   tabPanel("Decision Threshold",
        fluidRow(
          column(12,
            wellPanel(
              helpText("Introduction: Here is where you do calulations to determine the sample required to detect the 
                       differnece of interest between treatment arms. In other words, you're answering
                       the question, 'What difference do we need to see to scale this trial?'"),
              helpText("You can do this as an aboslute
                       (magnitude) change or as a percentage change.  The app will automatically show you
                       sample sizes for differences greater and less than the desired differnece to create
                       a menu of options for decision making."),
              checkboxInput(inputId = "percentage", label = "Percentage change?", value = FALSE),
              checkboxInput(inputId = "clustered", label = "Clustered design?", value = FALSE)
            )),
          column(4,
                 wellPanel(
              conditionalPanel(
                condition = "input.percentage == false & input.clustered == false",
                numericInput("sc_mean",
                             "Outcome average",
                             value = 100),
                numericInput("sc_stdev",
                             "Outcome standard deviation",
                             value = 50),
                numericInput("sc_diff_m",
                             "Absolute change over control",
                             value = 10)
              ), #conditional panel1
              conditionalPanel(
                condition = "input.percentage == false & input.clustered == true",
                numericInput("sc_mean_c",
                             "Outcome average",
                             value = 100),
                numericInput("sc_stdev_c",
                             "Outcome standard deviation",
                             value = 50),
                numericInput("sc_diff_mc",
                             "Absolute change over control",
                             value = 10),
                sliderInput("ICC_c", 
                            "Intra-cluster Correlation", 
                            value = 0.1, min = 0, max = 1),
                numericInput("sc_clustN_c",
                             "Average cluster size",
                             value = 10)
              ),
              conditionalPanel(
                condition = "input.percentage == true & input.clustered == false",
                numericInput("sc_mean_p",
                             "Outcome average",
                             value = 100),
                numericInput("sc_stdev_p",
                             "Outcome standard deviation",
                             value = 50),
                numericInput("sc_diff_p",
                             "Percentage (%) change over control",
                             value = 10)
              ),
              conditionalPanel(
                condition = "input.percentage == true & input.clustered == true",
                numericInput("sc_mean_pc",
                             "Outcome average",
                             value = 100),
                numericInput("sc_stdev_pc",
                             "Outcome standard deviation",
                             value = 50),
                numericInput("sc_diff_pc",
                             "Percentage (%) change over control",
                             value = 10),
                sliderInput("ICC_pc", 
                            "Intra-cluster Correlation", 
                            value = 0.1, min = 0, max = 1),
                numericInput("sc_clustN_pc",
                             "Average cluster size",
                             value = 10)
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
                       helpText("Introduction: Power calculations for when you have a set sample size and 
                                  you're trying to understand the type of effect you're 
                                likely to detect with your trial. This reports changes
                                in terms of standardized effect sizes. This app will automatically show
                                results for larger and smaller sample sizes to facilitate
                                decision making."),
                       helpText("If you do not enter a standard deviation, this page will only show you the Cohen's d
                                effect size. If you also enter the standard deviation of the outcome, the effect size 
                                will be converted into an average difference between treatment and control.")
                       )),
              column(4,
                     wellPanel(
                         numericInput("mde_n",
                                      "Available sample size per treatment arm",
                                      value = 100),
                         numericInput("mde_sd",
                                      "Standard deviation of outcome",
                                      value = NULL),
                         downloadButton("downloadMde", "Download Result")
                     ) #wellpanel layout
              ), # column close
              column(8,
                     wellPanel(
                       plotOutput("mde_plot"),
                       dataTableOutput("mde_table")
                     )
                  )
              )
        ),# tab panel close
    # new tabpanel of the PPV work?
  tabPanel("Positive Predictive Value",
           fluidRow(
             column(12,
                    wellPanel(
                      helpText("The purpose of this tab is place p-values in more context.
                               P-values help us understand the likelihood of seeing an effect as
                               great or greater than the observed effect but they don't tell us about
                               the overall likelihood of any difference being present. Positive predictive
                               value places observed p-values in the context of how likely it was to 
                               see an effect in the first place. The pre trial (a priori) expectation helps
                               adjust the the p-value to reflec the true false positive rate.")
                      )),
             column(4,
                    wellPanel(
                      sliderInput("priorOdds", label=h5("On a scale from 0 to 1, how likely is the hypothesis to be true?"),
                                  min = 0, max = 1, value = 0.25, step = 0.05),
                      sliderInput("power", label=h5("Select the power of the test"),
                                  min = 0, max = 1, value = 0.8, step = 0.01),
                      sliderInput("alpha", label=h5("Select observed alpha level"),
                                  min = 0, max = 1, value = 0.05, step = 0.01),
                      # selectInput("alpha", label = h5("Select alpha"),
                      #             choices = c("0.1", "0.05", "0.01"), selected = "0.05"),
                      br(),
                      actionButton("makeGraph", "Generate graph")
                    ) #wellpanel layout
             ), # column close
             column(8,
                    wellPanel(
                      plotOutput("ppv_plot"),
                      br(),
                      htmlOutput("ppv_percent")
                             )
                    )
                ) #fluidRow
        ),# tab panel close      
  tabPanel("Alpha Power Matrix",
           fluidRow(
             column(12,
                    wellPanel(
                      helpText("The purpose of this tab is to present the menu of feasible alpha
                               and power combinations available to trial designers given a difference
                               of interest and certain budget (sample size). This tab is most useful if trial designers
                               come with a clear picture of what an appropriate budget and sample size are
                               given the question of interest.")
                      )),
             column(4,
                    wellPanel(
                      numericInput("ap_mean",
                                   "Outcome average",
                                   value = 150),
                      numericInput("ap_stdev",
                                   "Outcome standard deviation",
                                   value = 80),
                      numericInput("ap_diff",
                                   "Absolute change over control",
                                   value = 25),
                      numericInput("ap_samp",
                                   "Available sample size",
                                   value=300),
                      br(),
                      actionButton("makeAlphaPower", "Generate graph")
                    ) #wellpanel layout
             ), # column close
             column(8,
                    wellPanel(
                      plotOutput("ap_plot", click = "plot_click"),
                      br(),
                      htmlOutput("ap_text")
                    )
             )
                      ) #fluidRow
             )# tab panel close      
)



server <- function(input, output, session) {
  
  # for decision thresholds
  # the same
  changes = sort(c(seq(-10,10, by=2), 1))
  
  # need to keep this simple 
  # this doesn't need to adjust for clustered designs because regardless I'm 
  # just feeding this into the function to look at different differences
  differs <- reactive({
    if(input$percentage==FALSE & input$clustered==FALSE){
      req(input$sc_mean, input$sc_stdev, input$sc_diff_m)
      toTest = input$sc_diff_m + changes
      toTest = toTest[toTest>0]
      return(toTest)
    } else if(input$percentage==FALSE & input$clustered==TRUE){
      req(input$sc_mean_c, input$sc_stdev_c, input$sc_diff_mc)
      toTest = input$sc_diff_mc + changes
      toTest = toTest[toTest>0]
      return(toTest)
    } else if(input$percentage==TRUE & input$clustered==FALSE){
      req(input$sc_mean, input$sc_stdev, input$sc_diff_p)
      toTest = ((input$sc_diff_p + changes)/100) * input$sc_mean_p
      toTest = toTest[toTest>0]
      return(toTest)
    } else if(input$percentage==TRUE & input$clustered==TRUE){
      req(input$sc_mean_pc, input$sc_stdev_pc, input$sc_diff_pc)
      toTest = ((input$sc_diff_pc + changes)/100) * input$sc_mean_pc
      toTest = toTest[toTest>0]
      return(toTest)
    }
})

  
  # the actual power calculations
  sampTab <- reactive({
    if(input$percentage==FALSE & input$clustered == FALSE){
      req(input$sc_mean, input$sc_stdev)
      dat = round(dat_MDE(input$sc_mean, input$sc_stdev, differs()),0)
    } else if(input$percentage==FALSE & input$clustered==TRUE){
      req(input$sc_mean_c, input$sc_stdev_c, input$sc_clustN_c, input$ICC_c)
      dat = round(dat_MDE_clus(input$sc_mean_c, input$sc_stdev_c, input$sc_clustN_c,
                               input$ICC_c, differs()),0)
    } else if(input$percentage == TRUE & input$clustered==FALSE){
      req(input$sc_mean_p, input$sc_stdev_p)
      dat = round(dat_MDE(input$sc_mean_p, input$sc_stdev_p, differs()),0)
    } else if(input$percentage == TRUE & input$clustered == TRUE){
      req(input$sc_mean_pc, input$sc_stdev_pc, input$sc_clustN_pc, input$ICC_pc)
      dat = round(dat_MDE_clus(input$sc_mean_pc, input$sc_stdev_pc, input$sc_clustN_pc,
                               input$ICC_pc, differs()),0)
    }
  })
  
  output$sc_plot <- renderPlot({
    
    dat = sampTab()
    xval = differs()
    
    gph <- ggplot() +
        geom_point(aes(x = xval, y = dat[,1]), size=3, color="green", shape=1) +
        geom_line(aes(x = xval, y = dat[,1]), size=1.2, color="green") + 
        geom_point(aes(x = xval, y = dat[,2]), size=3, color="blue", shape=1) +
        geom_line(aes(x = xval, y = dat[,2]), size=1.2, color="blue") +
        geom_point(aes(x = xval, y = dat[,3]), size=3, color="red", shape=1) +
        geom_line(aes(x = xval, y = dat[,3]), size=1.2, color="red") +
        labs(title = "Differences between treatments vs. sample size", y = "Sample size",
                  x = "Differences between treatments") + 
        theme(plot.title = element_text(hjust = 0.5,size = 20),
              plot.subtitle = element_text(hjust=0.5)) +
        scale_x_continuous(limits = c(min(xval), max(xval)))
    
    return(gph)
   
  })
  
  output$sc_table <- renderDataTable({
    
    farmerCap = 'Values are the number of farmers per treatment arm'
    clustCap = 'Values are the number of clusters per treatment arm'
    
    if(input$percentage==FALSE & input$clustered == FALSE){
    dat = format(sampTab(), big.mark = ",")
    # changes to evaluate in the data
    dat = as.data.frame(cbind(differs(), dat))
    names(dat) = c("Differences", "Alpha is 0.01", "Alpha is 0.05", "Alpha is 0.1")
    
    datatable(dat, rownames = FALSE, 
              options = list(pageLength = nrow(dat)),
              selection = list(mode = "single", 
                               selected=which(dat$Differences == input$sc_diff_m),
                               target="row"),
              caption = farmerCap)
    
    } else if(input$percentage==FALSE & input$clustered == TRUE){
      dat = format(sampTab(), big.mark = ",")
      # changes to evaluate in the data
      dat = as.data.frame(cbind(differs(), dat))
      names(dat) = c("Differences", "Alpha is 0.01", "Alpha is 0.05", "Alpha is 0.1")
      
      datatable(dat, rownames = FALSE,  
                options = list(pageLength = nrow(dat)),
                selection = list(mode = "single", 
                                 selected=which(dat$Differences == input$sc_diff_mc),
                                 target="row"),
                caption = clustCap)
    } else if(input$percentage==TRUE & input$clustered==FALSE){
      dat = format(sampTab(), big.mark = ",")
      # changes to evaluate in the data >> show both magnitude and percentage
      percentage = ((input$sc_diff_p + changes)/100)[((input$sc_diff_p + changes)/100)>0]
      
      
      dat = as.data.frame(cbind(differs(), percentage, dat))
      names(dat) = c("Differences", "Percent Change", "Alpha is 0.01", "Alpha is 0.05", "Alpha is 0.1")
      
      datatable(dat, rownames = FALSE,  
                options = list(pageLength = nrow(dat)),
                selection = list(mode = "single", 
                                 selected=which(dat$Differences == ((input$sc_diff_p/100)*input$sc_mean_p)),
                                 target="row"),
                caption = farmerCap)
    } else if(input$percentage==TRUE & input$clustered==TRUE){
      dat = format(sampTab(), big.mark = ",")
      # changes to evaluate in the data
      percentage = ((input$sc_diff_pc + changes)/100)[((input$sc_diff_pc + changes)/100)>0]
      
      dat = as.data.frame(cbind(differs(), percentage, dat))
      names(dat) = c("Differences", "Percent Change", "Alpha is 0.01", "Alpha is 0.05", "Alpha is 0.1")
      
      datatable(dat, rownames = FALSE,  
                options = list(pageLength = nrow(dat)),
                selection = list(mode = "single", 
                                 selected=which(dat$Differences == ((input$sc_diff_pc)/100*input$sc_mean_pc)),
                                 target="row"),
                caption = clustCap)
    }
  })
  
  
  # download data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sampleSize ", Sys.Date(), ".csv", sep = "")
    },
    
      content = function(file) {
        write.csv(sampTab(), file, row.names = FALSE)
      }
  )
  
  output$nrequired <- renderUI({
    
    # i'm doing this so it's easier to update messages.
    normalMsg = c("To achieve 80%% power and an alpha of 0.05 for a decision threshold of %d you'll need %s farmers per treatment arm")
    clusterMsg = c("To achieve 80%% power and an alpha of 0.05 for a decision threshold of %d, you'll need %s, 
                       clusters per treatment arm and %s farmers total per treatment arm")
    percentMsg = c("To achieve 80%% power and an alpha of 0.05 for a decision threshold of %d%% you'll need 
                   %s farmers per treatment arm")
    bothMsg = c("To achieve 80%% power and an alpha of 0.05 for a decision threshold of %d%% you'll need 
                %s clusters per treatment arm and %s farmers total per treatment arm")
    

    if(input$percentage==FALSE & input$clustered == FALSE){
      dat = sampTab()
      target = input$sc_diff_m
      samp_out = format(round(dat[differs()==input$sc_diff_m,2],0), big.mark = ",")
      
      str1 = sprintf(normalMsg, target, samp_out)
      
    } else if(input$percentage == FALSE & input$clustered == TRUE) {
      dat = sampTab()
      target = input$sc_diff_mc
      
      samp_out = round(dat[differs()==input$sc_diff_mc,2],0)
      farmersArm = format(samp_out * input$sc_clustN_c, big.mark = ",")
      samp_out = format(samp_out, big.mark = ",")
      
      str1 = sprintf(clusterMsg, target, samp_out, farmersArm)
      
    } else if(input$percentage == TRUE & input$clustered == FALSE) {
      dat = sampTab()
      target = input$sc_diff_p
      
      samp_out = format(round(dat[differs()==((input$sc_diff_p/100)*input$sc_mean_p),2],0), big.mark = ",")
      
      str1 = sprintf(percentMsg, target, samp_out)
      
    } else if(input$percentage == TRUE & input$clustered == TRUE) {
      dat = sampTab()
      target = input$sc_diff_pc
      
      samp_out = round(dat[differs()==((input$sc_diff_pc)/100*input$sc_mean_pc),2],0)
      farmersArm = format(samp_out * input$sc_clustN_pc, big.mark = ",")
      samp_out = format(samp_out, big.mark = ",")
      
      str1 = sprintf(bothMsg, target, samp_out, farmersArm)
      
    }
  })
  
  # minimum detectable effect section of server
  # add in look at getting mean if standard deviation is provided
  mde_changes = seq(0.1, 2, by=0.1)
  
  nOps = reactive({
    req(input$mde_n)
    nOps = input$mde_n * mde_changes
  })
  
  mde_dat <- reactive({
    req(input$mde_n)
    dat = do.call(rbind, lapply(nOps(), function(x){
      res = posthoc_mde(x)
      return(res)
    }))
    
    dat = data.frame(dat, n = nOps()) %>%
      setNames(c("a1", "a5", "a10", "n"))
    return(dat)
  })
  
  
  
  # minimum detectable effect
  output$mde_plot <- renderPlot({
    req(input$mde_n)
    
    #dat = data.frame(d = unlist(lapply(nOps(), posthoc_mde)), n = nOps())
    
    ggplot() + 
      geom_line(data = mde_dat(), aes(x = n, y = a1), size = 1.2, color="green") +
      geom_point(data=mde_dat()[which(mde_changes==1),], aes(x = n, y = a1), size=3) +
      geom_line(data = mde_dat(), aes(x = n, y = a5), size = 1.2, color="blue") +
      geom_point(data=mde_dat()[which(mde_changes==1),], aes(x = n, y = a5), size=3) +
      geom_line(data = mde_dat(), aes(x = n, y = a10), size = 1.2, color="red") +
      geom_point(data=mde_dat()[which(mde_changes==1),], aes(x = n, y = a10), size=3) +
      labs(x = "Sample size", y = "Minimum detectable effect (d)", 
           title = "Miniumum detectable effect for given sample size") + 
      theme(plot.title = element_text(hjust = 0.5, size = 20))
      
  })
  
  # do all transformations separate from datatable
  mde_tab <- reactive({
    if(is.null(input$mde_sd)){
      
      dat = melt(mde_dat(), id.vars = "n") %>%
        dplyr::mutate(variable = ifelse(variable == "a1", "0.1",
                                        ifelse(variable == "a5", "0.05", "0.01")),
                      value = round(value, 3))  
      
      names(dat) = c("Sample Size", "Alpha level less than", "Effect Size")
      dat = as.data.frame(dat)
    
    } else {
      
      dat = melt(mde_dat(), id.vars = "n") %>%
        dplyr::mutate(variable = ifelse(variable == "a1", "0.1",
                                        ifelse(variable == "a5", "0.05", "0.01")),
                      value = round(value, 3))
      
      dat$mu = dat$value * input$mde_sd
      names(dat) = c("Sample Size", "Alpha level less than", "Effect Size", "Average Difference")
      dat = as.data.frame(dat)
    }
  })
  
  # mde output table
  output$mde_table <- renderDataTable({
    
    datatable(mde_tab(), rownames = FALSE, selection = list(mode = "single", 
                                                      selected=which(mde_changes==1),
                                                      target="row"),
              filter = 'top',
              options = list(pageLength = 20))
   
  })
  
  # download mde result
  output$downloadMde <- downloadHandler(
    filename = function() {
      paste("mdeResult ", Sys.Date(), ".csv", sep = "")
    },
    
    content = function(file) {
      write.csv(mde_tab(), file, row.names = FALSE)
    }
  )
  
  # ppv
  plotTiles <- eventReactive(input$makeGraph, {
    
    # alpha <- ifelse(input$alpha=="0.1", 0.1, 
    #                 ifelse(input$alpha=="0.05", 0.05, 0.01))
    
    trueHypo <- round(100 * input$priorOdds,0)
    falseHypo <- round(100 - trueHypo,0)
    
    # falses
    falsePos <- round(falseHypo * input$alpha,0)
    falseHypo <- round(falseHypo - falsePos,0)
    
    # truths
    falseNeg <- round(trueHypo * (1 - input$power),0)
    truePos <- round(trueHypo * input$power,0)
    
    nrows <- 10
    df <- expand.grid(y = 1:nrows, x = 1:nrows)
    
    var = c(rep("1. False Hypotheses", falseHypo), rep("2. False Positives", falsePos),
            rep("3. False Negatives", falseNeg), rep("4. True Positives", truePos))
    categ_table <- round(table(var) * ((nrows*nrows)/(length(var))))
    categ_table
    
    df$category <- factor(rep(names(categ_table), categ_table))
    
    ggplot(df, aes(x = x, y = y, fill = category)) + 
      geom_tile(color = "black", size = 0.5) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
      scale_fill_brewer(palette = "Set3") +
      #coord_equal() +
      labs(title="Hypotheses: Why power matters for finding true effects", 
           subtitle="Being powerful is good",
           caption="Source: the Economist",
           x = "", y = "") + 
      theme(legend.position = "bottom")
    
  })
  
  trueAlpha <- eventReactive(input$makeGraph, {
    
    # alpha <- ifelse(input$alpha=="0.1", 0.1, 
    #                 ifelse(input$alpha=="0.05", 0.05, 0.01))
    
    trueHypo <- round(100 * input$priorOdds,0)
    falseHypo <- round(100 - trueHypo,0)
    
    # falses
    falsePos <- round(falseHypo * input$alpha,0)
    falseHypo <- round(falseHypo - falsePos,0)
    
    # truths
    falseNeg <- round(trueHypo * (1 - input$power),0)
    truePos <- round(trueHypo * input$power,0)
    
    percentage <- paste(round(falsePos / (truePos + falsePos) * 100,1), "%")
    res <- paste0("True false positive rate: ", percentage)
    
  })
  
  # and ouput those results
  output$ppv_plot <- renderPlot({
    plotTiles()
  })
  
  output$ppv_percent <- renderText({
    trueAlpha()
  })
  
  
  
  apDat <- eventReactive(input$makeAlphaPower, {
    samp1 <- rnorm(n=1000, mean = input$ap_mean, sd=input$ap_stdev)
    #this is a better version if you can understand it:
    samp2 <- samp1 + rnorm(length(samp1), input$ap_diff, input$ap_diff/10) #add some noise
    inp <- abs(cohen_d(samp1, samp2)$effectsi)
    
    pvals <- seq(0.05, 0.5, by=0.05)
    powerVals <- seq(0.5, 0.95, by=0.05)
    
    res <- do.call(rbind, lapply(pvals, function(pval){
      
      do.call(rbind, lapply(powerVals, function(powerNum){
        
        result = tryCatch({
          n <- pwr.2p.test(h=inp, sig.level=pval, power = powerNum, n=NULL)$n
        }, warning = function(w) {
          n <- NA
        }, error = function(e){
          n <- NA
        })
        
        figs = data.frame(pval, powerNum, result)
        return(figs)
        
      }))
    }))
    
    #make graph
    nrows <- length(pvals)
    df <- expand.grid(y = 1:nrows, x = 1:nrows)
    
    # var = c(rep("1. False Hypotheses", falseHypo), rep("2. False Positives", falsePos),
    #         rep("3. False Negatives", falseNeg), rep("4. True Positives", truePos))
    # categ_table <- round(table(var) * ((nrows*nrows)/(length(var))))
    # categ_table
    df <- cbind(df, res)
    df$result <- round(df$result)
    df$check <- ifelse(df$result <= input$ap_samp, TRUE, FALSE) # TRUE IS GOOD, FALSE IS BAD
    
    #df$category <- factor(rep(names(categ_table), categ_table))
    df$category <- factor(ifelse(df$check==TRUE, "Current sample sufficient", 
                                 ifelse(df$check==FALSE, "Need bigger sample", "NA")))
    return(df)
  })
  
  # alpha power matrix
  plotAlphaPower <- eventReactive(input$makeAlphaPower, {
    
    df <- apDat()
    
    pvals <- seq(0.05, 0.5, by=0.05)
    powerVals <- seq(0.5, 0.95, by=0.05)
    
    ggplot(df, aes(x = x, y = y, fill = category)) + 
      geom_tile(color = "black", size = 0.5) +
      geom_text(label = round(df$res, 1)) +
      scale_x_continuous(expand = c(0, 0), breaks = unique(df$x), labels = pvals) +
      scale_y_continuous(expand = c(0, 0), trans = 'reverse', breaks = unique(df$y), labels = powerVals) +
      scale_fill_brewer(palette = "Set3") +
      #coord_equal() +
      labs(title="Trade-offs of alpha and power", 
           caption="Source: our brainz",
           x = "p-values (decreasing)", y = "Power", category = "Labels") + 
      theme(legend.position = "bottom")
    
    
    
  })
  
  
  #alpha power plot
  output$ap_plot <- renderPlot({
    plotAlphaPower()
  })
  
  # for clicking on the plot
  clickDat <- reactive({
    clickCoords = data.frame(x = round(as.numeric(input$plot_click$x)), y=round(as.numeric(input$plot_click$y)))
    df = apDat()
    val = df[df$x==clickCoords$x & df$y==clickCoords$y, c("result", "pval", "powerNum")]
    return(val)
  })
  
  #alpha power text
  output$ap_text <- renderUI({

    dat = clickDat()
    dat$pval = round(dat$pval,2)
    #dat$wrong = dat$pval * (1-dat$powerNum)*100
    #normalWrong = (0.05*0.2)*100
    
    #dat$increase = (dat$wrong - normalWrong)/normalWrong
    boiler = "This combination is an alpha of %s and power of %d%% requires a sample of %d. This means we have a %s%% chance of a false positive and a %s%% chance of a false negative."
    str = sprintf(boiler, as.character(dat$pval), dat$powerNum*100, dat$result, as.character(dat$pval*100), as.character((1-dat$powerNum)*100))
    return(str)

  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)


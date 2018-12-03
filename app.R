library(ggplot2)
library(magrittr)
library(rstudioapi)
library(KMsurv)
library(survival)
library(survminer)
library(shiny)
library(ggplot2)
library(plyr)
library(purrr)
library(taRifx)
library(plotrix)
library(tidyverse)

tx <- readRDS(file = "shiny_data.rds")

## Functions

# This Functions converts P-Values of less than 0.001 to show "<0.001" as to not display 0
print_p <- function(a_vector) {
  ifelse(a_vector <= 0.001, "<0.001", as.character(a_vector))
}

# I don't actually think this is being used in this app
p_gsub <- function(a_vec, to_repl, repl_with) {
  stopifnot(length(to_repl) == length(repl_with))
  reduce2(to_repl, 
         repl_with, 
         function(a, tr, rw) gsub(tr, rw, a), 
         .init = a_vec
         )
}

# Removes rows for missing data in selected columns, used for Survival analysis
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Calculates Mean and Standard Deviation of Variable by Category
data_mean_by_cat <- function(data, varname, groupnames){
  require(plyr)
  require(plotrix)
  summary_func <- function(x, col){
    c(Mean = mean(x[[col]], na.rm=TRUE),
      SD = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}

# Amended dataset for survival analysis removing patients with missing Vital Status and Follow-Up Information
tx_surv <- completeFun(tx, c("PUF_VITAL_STATUS","DX_LASTCONTACT_DEATH_MONTHS"))

# Categorical Predictors
predictors <- 
  c(`Facility Type`= "Facility_Type", 
    "Sex",
    "Race",
    "Ethnicity",
    `Charlson/Deyo Score` = "Charlson_Deyo_Score",
    `Facility Location`= "Facility_Location",
    `Insurance Status` = "Insurance_Status" ,
    `Urban/Rural Status` = "Urban_Rural_Status",
    `Tumor Behavior` = "Behavior",
    `AJCC Clinical Stage Group` = "AJCC_Clinical_Stage_Group",
    `AJCC Pathologic Stage Group` = "AJCC_Pathologic_Stage_Group",
    `Bone Metastases` = "Bone_Metastases",
    `Lung Metastases` = "Lung_Metastases",
    `Radiation Therapy` = "Radiation_Therapy",
    "Chemotherapy",
    `30 Day Mortality` = "Thirty_Day_Mortality",
    `90 Day Mortality` = "Ninety_Day_Mortality"
    )

# Continuous Predictors
predictors.cont <- c("Age",
                    `Surgical Inpatient Stay` = "Surgical_Inpatient_Stay_Days",
                    `Days Treatment Started from Diagnosis` = "Days_Treatment_Started_from_Diagnosis",
                    `Regional Radiation Dose` = "Regional_Radiation_Dose",
                    `Tumor Size` = "Tumor_Size")

# Javascript Array for Proportions Button
the_condition <- 
  "['Facility_Type', 
    'Sex',
    'Race',
    'Ethnicity',
    'Charlson_Deyo_Score',
    'Facility_Location',
    'Insurance_Status',
    'Urban_Rural_Status',
    'Behavior',
    'AJCC_Clinical_Stage_Group',
    'AJCC_Pathologic_Stage_Group',
    'Bone_Metastases',
    'Lung_Metastases',
    'Radiation_Therapy',
    'Chemotherapy',
    'Thirty_Day_Mortality',
    'Ninety_Day_Mortality'].includes(input.predictor)"

# UI
ui <- fluidPage(
       
  
  
        navbarPage("NCDB",
                   
                   tabPanel("Descriptive Statistics",    
                            sidebarLayout(
                              sidebarPanel(
                                
                                selectInput('location', 'Location',
                                            c(All = 'All', unique(tx$Cancer_site))),
                                selectInput('predictor', 'Variable of Interest', 
                                            c(predictors, predictors.cont)
                                            ),
                                selectInput('facet', "Alt. sub-plots", 
                                            c(None = ".", 
                                              predictors)),
                                conditionalPanel(
                                  condition = the_condition,
                                  checkboxInput('counts', 'Proportions instead of counts')
                                  )
                              ),
                              mainPanel(
                                h4(tableOutput("desc_title")),h1(plotOutput("plot", height = "550px"))
                                ,
                                fluidRow(
                                  column(width = 12, tableOutput("desc"))
                                
                                )))),
                                
                              
                   
                  
                   tabPanel("Survival Analysis - Kaplan-Meier Survival Graph",
         sidebarLayout(
           sidebarPanel(
             h3("Survival Graph"),
             
             selectInput('site', 'Location', 
                         c(None = 'All', unique(tx_surv$Cancer_site))),
                         
             selectInput('sur_var', 'Factor of Survival', 
                         c(predictors)),
             
             sliderInput("xvalue", 'Survival Months = ',
                         value=100, min=1, max=max(tx_surv$DX_LASTCONTACT_DEATH_MONTHS))
           )
           ,
           mainPanel(
           h3(textOutput("caption")),plotOutput("plot1")
           , tableOutput("center"), tableOutput("left")
           
           
           )))))



server <- function(input, output, session) {
  
  dataset <- reactive({
    
    all_true <- rep_len(TRUE, nrow(tx))
    if (input$location != 'All') {
      keep_these_sites <- (tx[["Cancer_site"]] == input$location)
    } else
      keep_these_sites <- all_true
   
     tx[keep_these_sites,]
    
  })
  
  counts.df <- reactive({
    
    if ((input$predictor %in% c(predictors)) & input$facet == '.'){
      counts <- plyr::count(dataset(), c(input$predictor))
      counts$label <- ((counts$freq/(sum(counts$freq)) * 100))
      counts$percent = paste0(sprintf("%.0f", counts$label), "%")
      
    } 
    
    if ((input$predictor %in% c(predictors)) & input$facet != '.'){
      counts <- plyr::count(dataset(), c(input$predictor, c(input$facet)))
      counts <- ddply(counts, .(counts[,1]), transform, label = freq/sum(freq) * 100)
      counts$percent = paste0(sprintf("%.0f", counts$label), "%")
      
      pre.1 <- dataset()[, c(input$predictor)]
      fac.1 <- dataset()[, c(input$facet)]
      dtb <- cbind(pre.1,fac.1)
      tbl <- table(dtb)
      pval <- chisq.test(tbl)
      pval$p.value <- round(pval$p.value, digits = 3)
      pval$p.value <- print_p( pval$p.value)
      p_table <- as.data.frame(pval$p.value)
      names(p_table) <- c("P Value")
      p_table <- as.data.frame(p_table)
      
      counts <- cbind(counts, p_table, row.names = NULL)
      counts[2:nrow(counts),"P Value"] <- NA
      } 
    
    if ((input$predictor %in% c(predictors.cont)) & input$facet != '.'){
      counts <- data_mean_by_cat(dataset(), c(input$predictor), c(input$facet))
      
      out.1 <- dataset()[, c(input$predictor)]
      out.2 <- as.numeric(unlist(out.1))
      fac.1 <- dataset()[, c(input$facet)]
      fac.2 <- as.factor(unlist(fac.1))
      
      fit <- summary(aov(out.2 ~ fac.2))
      sum_test <- unlist(fit)
      names(sum_test)
      p.value <- sum_test["Pr(>F)1"]
      p.value <- round(p.value, digits = 3)
      p.value <- print_p(p.value)
      p_table <- as.data.frame(as.matrix(p.value))
      colnames(p_table)[colnames(p_table)=="V1"] <- "P Value"
      
      
      counts <- cbind(counts, p_table, row.names = NULL)
      counts[2:nrow(counts),"P Value"] <- NA
      
    }
    
    if ((input$predictor %in% c(predictors.cont)) & input$facet == '.'){
      
      a <- dataset()[, c(input$predictor)]
      x <- as.numeric(unlist(a))
      counts <- c(
                  mean = mean(x, na.rm = TRUE) 
                 , sd = sd(x, na.rm = TRUE))
      counts <- matrix(counts, ncol=2, byrow=TRUE)
      counts <- as.data.frame(counts, stringsAsFactors=FALSE)
      names(counts) <- c("Mean", "SD")
      
  
      
    }
    counts
    
  })
  
  output$plot <- renderPlot({
      
      if ((input$predictor %in% c(predictors)))
        
        p <- ggplot(counts.df(), 
             aes_string(x = input$predictor, y = "freq" )
      ) + 
      geom_bar(stat="identity",colour="black", position = position_dodge())
    
      if ((input$predictor %in% c(predictors.cont)))
        p <- 
          ggplot(dataset(), 
                 aes_string(x = input$predictor)
          ) + 
          geom_histogram(position = position_dodge(), binwidth = 1) 
      
      if (input$counts & (input$predictor %in% c(predictors)) )
        p <- ggplot(counts.df(), 
                    aes_string(x = input$predictor, y = "label" )
        ) + 
          geom_bar(stat="identity",colour="black", position = position_dodge()) 
      
      if ((input$predictor %in% c(predictors)) & input$facet != '.')
      p <- p + aes_string(fill = input$facet)
      
      if ((input$predictor %in% c(predictors.cont)) & input$facet != '.')
        p <- p + facet_wrap(input$facet, scales = "free_y")
      
    p + theme_bw()
    
  })
  
  output$desc_title <- renderText({
    
    
    if ((input$predictor %in% c(predictors)) & input$facet != '.')
      desc_title <- paste("Distribution of", input$predictor,"by",input$facet, sep = "\n")
    
    if ((input$predictor %in% c(predictors)) & input$facet == '.')   
      desc_title <- paste("Distribution of", input$predictor, sep = "\n")

    if ((input$predictor %in% c(predictors.cont)) & input$facet == '.')
      desc_title <- paste("Distribution of", input$predictor,sep = "\n")
        
    if ((input$predictor %in% c(predictors.cont)) & input$facet != '.')
      desc_title <- paste("Distribution of", input$predictor,"by",input$facet, sep = "\n")
    
    desc_title
    
  })
  
  output$desc <- renderTable({
    
    if ((input$predictor %in% c(predictors)) & input$facet != '.'){
      

    df <- counts.df()  %>% mutate(zero_p = paste(freq," ","(",percent,")",sep=""))
    keep <- c(input$predictor,input$facet,"zero_p")
    df.l <- df[keep]
    
    library(reshape2)
    # df.z<- recast(df.l, as.formula(paste(paste(input$facet)," + zero_p  ~ ",paste(input$predictor))))
    # 
    # df1 <- subset(df.z, select = -c(zero_p))
    
    final <-  spread(df.l, paste(input$predictor), zero_p)
    
    library(data.table)
    # final <- setDT(df1)[, lapply(.SD, na.omit), by = df1[,1]]
    leng <- nrow(final)
    p.val <- df$`P Value`
    length(p.val) <- leng
    desc_table <- cbind(final,p.val)
    }
    
    if ((input$predictor %in% c(predictors)) & input$facet == '.')   
    desc_table <- as.data.frame(counts.df() [c(input$predictor, "freq", "percent")])
    
    if ((input$predictor %in% c(predictors.cont)))
      desc_table <- as.data.frame(counts.df())
      
    

    desc_table
  } , na = "")

  
  selectedData <- reactive({
    
    all_true.1 <- rep_len(TRUE, nrow(tx_surv))
    
    if (input$site != 'All') {
      keep_these_sites.1 <- (tx_surv[["Cancer_site"]] == input$site)
    } else
      keep_these_sites.1 <- all_true.1
    
    tx_surv[keep_these_sites.1,c(input$sur_var,'DX_LASTCONTACT_DEATH_MONTHS','PUF_VITAL_STATUS')]
  })
  
  output$caption <- renderText({
    paste("Survival Graph of", input$sur_var, sep = "\n")
  })
  
  runSur <- reactive({
    survfit(as.formula(paste("Surv(DX_LASTCONTACT_DEATH_MONTHS,PUF_VITAL_STATUS) ~ ",paste(input$sur_var))),data=selectedData())
  })
  
  output$plot1 <- renderPlot({
    
    plot(runSur(), 
         col=c("red","sky blue","green","purple","orange","yellow"), xlab="Months", ylab="S(t)")
    lLab <- gsub(pattern = paste(input$sur_var,"=",sep = ""), replacement = "",(names(runSur()$strata)))
    legend(
      "topright",
      cex=0.9,
      legend = lLab,
      levels(selectedData()),
      fill= c("red","sky blue","green","purple","orange","yellow"))
    abline(v=input$xvalue,col=1,lty=2)
  })
  
  output$center <- renderTable({
    tbl <- as.data.frame(summary(runSur(), times=input$xvalue) [c("strata","time","surv")])
    tbl[["strata"]] <- gsub(pattern = paste(input$sur_var,"=",sep = ""), replacement = "", x = tbl[["strata"]])
    colnames(tbl)[colnames(tbl)=="surv"] <- "Survival Probability"
    tbl
  }) 
  
  
  output$left <- renderTable({
    btm.tbl <-as.data.frame(surv_median(runSur()) [c("strata", "median", "upper", "lower")])
    btm.tbl[["strata"]] <- gsub(pattern = paste(input$sur_var,"=",sep = ""), replacement = "", x = btm.tbl[["strata"]])
    btm.tbl$`Median Survival (95% CI)` <-  paste(btm.tbl$median," ","(",btm.tbl$lower,"-",btm.tbl$upper,")",sep="")
    btm.tbl.f<- btm.tbl[c("strata","Median Survival (95% CI)")]
    btm.tbl.f
    
  })
  cox <- reactive({
    coxph(as.formula(paste("Surv(DX_LASTCONTACT_DEATH_MONTHS,PUF_VITAL_STATUS) ~ ",paste(input$sur_var))),data=selectedData())
  })
  
}



shinyApp(ui, server)



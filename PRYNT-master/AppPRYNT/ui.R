ui <- fluidPage(
  
  # App title ----
  mainPanel(h1("Application PRYNT"),h5("PRioritization bY protein NeTwork")),
  mainPanel(hr(),
            h4("Prioritize disease candidates from urinary protein profile"),br(),br() ),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    # sidebarPanel(
      verticalLayout(
        wellPanel(
          fluidRow(
            column(2),
            column(6,
              textAreaInput("caption", "List of deregulated proteins (one of each line)", "", width = "600px",height = "200px"),
              fluidRow(
                column(2,"Not found : "),
                column(6,
                verbatimTextOutput("value"),tags$head(tags$style("#value{color: red;font-size: 20px;font-style: bold;}")),
                )
              ),hr(),
              fluidRow(
                column(3),
                column(3,
                  actionButton("confirmdatabutton",h4("Run PRYNT"),align="center", width ='300px')
                )
              )
            )
          )
        )
    
    
),
    
    # Main panel for displaying outputs ----
    mainPanel(
      p(conditionalPanel(condition ="input.confirmdatabutton!=0" ,
            h2("Results"),
            # withSpinner(ui_element = 
            dataTableOutput('table') %>% withSpinner(color="#0dc5c1",type = 8),align="center"),

            # ,color="#0dc5c1",type = 8,))
            p(conditionalPanel(condition ="input.confirmdatabutton!=0" ,
                    downloadButton("downloadtable","Download"),align="center", width ='200px')
            )
            

      )
      # dataTableOutput("table")
    )
  )
)

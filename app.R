library(shiny)
library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)
library(splines)
library(Lahman)
library(rsconnect)
# install.packages(c("fastmap", "bslib", "htmltools", "shiny"))
# install.packages("htmltools")
# install.packages("shiny")

# ============================================================
# LOAD + PREP DATA
# For now, this includes your project code directly.
# Later, you can move this into a separate script or saved files.
# ============================================================

# ---------- Statcast data ----------
statcast_raw <- read.csv("C:/Users/Alicia/Downloads/SMGT490/stats.csv") %>%
  clean_names()


statcast_raw <- statcast_raw %>%
  mutate(
    pitcher_id   = player_id,
    pitcher_name = last_name_first_name,
    year         = as.integer(year),
    age          = as.integer(player_age),
    games        = p_starting_p + p_game_in_relief,
    start_share  = ifelse(games > 0, p_starting_p / games, NA_real_),
    relief_share = ifelse(games > 0, p_game_in_relief / games, NA_real_),
    era_type     = "Statcast"
  )

statcast_clean <- statcast_raw %>%
  filter(!is.na(p_era), !is.na(age), !is.na(pitcher_id)) %>%
  filter(p_formatted_ip >= 50) %>%
  filter(n >= 200) %>%
  filter(!is.na(fastball_avg_speed)) %>%
  mutate(
    v_z       = as.numeric(scale(fastball_avg_speed)),
    spin_z    = as.numeric(scale(ff_avg_spin)),
    pitches_z = as.numeric(scale(pitch_count)),
    ip_z      = as.numeric(scale(p_formatted_ip)),
    whiff_z   = as.numeric(scale(whiff_percent)),
    bb_z      = as.numeric(scale(bb_percent)),
    k_z       = as.numeric(scale(k_percent))
  )

pitcher_peaks <- statcast_clean %>%
  arrange(pitcher_id, p_era, desc(p_formatted_ip), desc(pitch_count)) %>%
  group_by(pitcher_id, pitcher_name) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    pitcher_name,
    pitcher_id,
    peak_year = year,
    peak_age  = age,
    peak_era  = p_era,
    ip        = p_formatted_ip,
    k_pct     = k_percent,
    bb_pct    = bb_percent,
    fb_v      = fastball_avg_speed
  )

m1_statcast <- lmer(
  p_era ~ ns(age, df = 4) +
    v_z + spin_z + pitches_z + ip_z + whiff_z + start_share +
    factor(year) +
    (1 | pitcher_id),
  data = statcast_clean,
  REML = FALSE
)

baseline_grid <- expand_grid(
  age         = 21:40,
  v_z         = 0,
  spin_z      = 0,
  pitches_z   = 0,
  ip_z        = 0,
  whiff_z     = 0,
  start_share = median(statcast_clean$start_share, na.rm = TRUE),
  year        = median(statcast_clean$year, na.rm = TRUE)
)

baseline_grid$pred_era <- predict(m1_statcast, newdata = baseline_grid, re.form = NA)

m2_velocity <- lmer(
  p_era ~ ns(age, df = 4) * v_z +
    spin_z + pitches_z + ip_z + whiff_z + start_share +
    factor(year) +
    (1 | pitcher_id),
  data = statcast_clean,
  REML = FALSE
)

velocity_grid <- expand_grid(
  age         = 21:40,
  v_z         = c(-1, 0, 1),
  spin_z      = 0,
  pitches_z   = 0,
  ip_z        = 0,
  whiff_z     = 0,
  start_share = median(statcast_clean$start_share, na.rm = TRUE),
  year        = median(statcast_clean$year, na.rm = TRUE)
)

velocity_grid$pred_era <- predict(m2_velocity, newdata = velocity_grid, re.form = NA)

vol_data <- statcast_clean %>%
  arrange(pitcher_id, year) %>%
  group_by(pitcher_id) %>%
  mutate(
    era_next  = lead(p_era),
    age_next  = lead(age),
    d_era     = era_next - p_era,
    abs_d_era = abs(d_era)
  ) %>%
  ungroup() %>%
  filter(!is.na(abs_d_era), !is.na(age_next)) %>%
  mutate(age_mid = age)

m_volatility <- lmer(
  abs_d_era ~ ns(age_mid, df = 4) +
    v_z + spin_z + pitches_z + ip_z + whiff_z + start_share +
    factor(year) +
    (1 | pitcher_id),
  data = vol_data,
  REML = FALSE
)

vol_grid <- expand_grid(
  age_mid     = 21:40,
  v_z         = 0,
  spin_z      = 0,
  pitches_z   = 0,
  ip_z        = 0,
  whiff_z     = 0,
  start_share = median(vol_data$start_share, na.rm = TRUE),
  year        = median(vol_data$year, na.rm = TRUE)
)

vol_grid$pred_vol <- predict(m_volatility, newdata = vol_grid, re.form = NA)

# ---------- Historical data ----------
lahman_compare <- Pitching %>%
  left_join(People %>% select(playerID, birthYear), by = "playerID") %>%
  mutate(
    pitcher_id = as.character(playerID),
    year       = yearID,
    age        = yearID - birthYear,
    ip         = IPouts / 3,
    ERA        = (ER * 9) / ip,
    k_rate     = SO / ip,
    bb_rate    = BB / ip,
    era_type   = "Pre-Statcast"
  ) %>%
  filter(year >= 1985, year <= 2000) %>%
  filter(ip >= 50, !is.na(age), !is.na(ERA), is.finite(ERA))

m_hist <- lmer(
  ERA ~ ns(age, df = 4) +
    k_rate + bb_rate + ip +
    factor(year) +
    (1 | pitcher_id),
  data = lahman_compare,
  REML = FALSE
)

hist_grid <- expand_grid(
  age     = 21:40,
  k_rate  = mean(lahman_compare$k_rate, na.rm = TRUE),
  bb_rate = mean(lahman_compare$bb_rate, na.rm = TRUE),
  ip      = mean(lahman_compare$ip, na.rm = TRUE),
  year    = median(lahman_compare$year, na.rm = TRUE)
)

hist_grid$pred_era <- predict(m_hist, newdata = hist_grid, re.form = NA)

modern_compare <- statcast_clean %>%
  mutate(
    pitcher_id = as.character(pitcher_id),
    ERA        = p_era,
    ip         = p_formatted_ip,
    k_rate     = k_percent / 100,
    bb_rate    = bb_percent / 100,
    era_type   = "Statcast"
  ) %>%
  select(pitcher_id, year, age, ERA, ip, k_rate, bb_rate, era_type)

combined_data <- bind_rows(
  modern_compare,
  lahman_compare %>% select(pitcher_id, year, age, ERA, ip, k_rate, bb_rate, era_type)
) %>%
  mutate(
    era_type  = factor(era_type),
    k_rate_z  = as.numeric(scale(k_rate)),
    bb_rate_z = as.numeric(scale(bb_rate)),
    ip_z      = as.numeric(scale(ip))
  )

combined_model <- lmer(
  ERA ~ ns(age, df = 4) * era_type +
    k_rate_z + bb_rate_z + ip_z +
    factor(year) +
    (1 | pitcher_id),
  data = combined_data,
  REML = FALSE
)

combined_grid <- expand_grid(
  age       = 24:38,
  era_type  = c("Pre-Statcast", "Statcast"),
  k_rate_z  = 0,
  bb_rate_z = 0,
  ip_z      = 0,
  year      = median(combined_data$year, na.rm = TRUE)
)

combined_grid$pred_era <- predict(combined_model, newdata = combined_grid, re.form = NA)

# ============================================================
# SHINY APP HELPERS
# ============================================================

# Create searchable player choices from players actually in the dataset
player_choices <- statcast_clean %>%
  distinct(pitcher_id, pitcher_name) %>%
  arrange(pitcher_name) %>%
  mutate(
    pitcher_id = as.character(pitcher_id),
    label = paste0(pitcher_name, " (", pitcher_id, ")")
  )

# Small helper function for MLB headshot URL
get_headshot_url <- function(player_id) {
  paste0("https://content.mlb.com/images/headshots/current/60x60/", player_id, ".png")
}


# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body {
        background-color: #f5f7fb;
        font-family: 'Arial', sans-serif;
      }

      .app-title {
        background: linear-gradient(135deg, #0b1f3a, #1f5f99);
        color: white;
        padding: 25px;
        border-radius: 0 0 18px 18px;
        margin-bottom: 25px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
      }

      .app-title h1 {
        margin: 0;
        font-weight: 700;
      }

      .app-title p {
        margin-top: 8px;
        font-size: 16px;
        opacity: 0.9;
      }

      .card {
        background: white;
        border-radius: 16px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.08);
      }

      .metric-card {
        background: white;
        border-radius: 14px;
        padding: 15px;
        text-align: center;
        box-shadow: 0 3px 10px rgba(0,0,0,0.08);
        margin-bottom: 15px;
      }

      .metric-value {
        font-size: 26px;
        font-weight: bold;
        color: #12355b;
      }

      .metric-label {
        font-size: 13px;
        color: #666;
      }

      .player-card {
        display: flex;
        align-items: center;
        gap: 18px;
      }

      .player-img {
        width: 75px;
        height: 75px;
        border-radius: 50%;
        border: 3px solid #12355b;
        background-color: #eaeaea;
      }

      .player-name {
        font-size: 24px;
        font-weight: 700;
        margin-bottom: 3px;
      }

      .player-sub {
        color: #666;
        font-size: 14px;
      }

      .section-title {
        font-weight: 700;
        color: #12355b;
        margin-bottom: 15px;
      }
    "))
  ),
  
  div(
    class = "app-title",
    h1("MLB Pitcher Aging Curves Dashboard"),
    p("Interactive dashboard for exploring Statcast-era aging curves, historical comparisons, volatility, and individual pitcher profiles.")
  ),
  
  tabsetPanel(
    
    # ========================================================
    # TAB 1: OVERVIEW DASHBOARD
    # ========================================================
    tabPanel(
      "Dashboard",
      br(),
      
      sidebarLayout(
        sidebarPanel(
          div(
            class = "card",
            h4("Choose a Visualization"),
            selectInput(
              "plot_type",
              "Graph:",
              choices = c(
                "Peak Age Distribution",
                "Statcast Aging Curve",
                "Velocity Aging Curves",
                "Volatility Curve",
                "Pre-Statcast Aging Curve",
                "Era Comparison"
              )
            )
          )
        ),
        
        mainPanel(
          div(
            class = "card",
            plotOutput("mainPlot", height = "450px")
          ),
          
          div(
            class = "card",
            h4("Interpretation"),
            htmlOutput("summaryText")
          )
        )
      )
    ),
    
    
    # ========================================================
    # TAB 2: PLAYER EXPLORER
    # ========================================================
    tabPanel(
      "Player Explorer",
      br(),
      
      sidebarLayout(
        sidebarPanel(
          div(
            class = "card",
            h4("Search Player"),
            selectizeInput(
              "selected_player",
              "Choose a pitcher:",
              choices = setNames(player_choices$pitcher_id, player_choices$label),
              selected = player_choices$pitcher_id[1],
              options = list(
                placeholder = "Type a player name...",
                maxOptions = 20
              )
            )
          )
        ),
        
        mainPanel(
          
          div(
            class = "card",
            uiOutput("playerProfile")
          ),
          
          fluidRow(
            column(3, uiOutput("playerYearsCard")),
            column(3, uiOutput("playerPeakAgeCard")),
            column(3, uiOutput("playerPeakEraCard")),
            column(3, uiOutput("playerVelocityCard"))
          ),
          
          div(
            class = "card",
            h4(class = "section-title", "Season-by-Season Performance"),
            plotOutput("playerSeasonPlot", height = "380px")
          ),
          
          div(
            class = "card",
            h4(class = "section-title", "Player Season Table"),
            tableOutput("playerSeasonTable")
          )
        )
      )
    ),
    
    
    # ========================================================
    # TAB 3: ABOUT
    # ========================================================
    tabPanel(
      "About Project",
      br(),
      
      div(
        class = "card",
        h2("About This Project"),
        p(
          "This dashboard explores how Major League Baseball pitchers age over time, 
      with a focus on whether modern pitchers in the Statcast era reach peak 
      performance earlier or decline faster than pitchers from earlier eras."
        ),
        p(
          "The project connects traditional pitcher aging research with modern 
      Statcast-era variables such as velocity, spin rate, pitch count, whiff rate, 
      innings pitched, and workload intensity."
        )
      ),
      
      div(
        class = "card",
        h3("Research Question"),
        p(
          strong("How have MLB pitcher aging patterns changed in the Statcast era, 
      and how do factors such as velocity, workload, and strikeout ability influence 
      pitcher peak performance and career decline?")
        )
      ),
      
      div(
        class = "card",
        h3("Introduction"),
        p(
          "Pitching in Major League Baseball has changed dramatically with the rise of 
      Statcast tracking technology. Since 2015, teams have been able to measure 
      pitching performance in much greater detail, including velocity, spin rate, 
      release point, pitch movement, and workload patterns."
        ),
        p(
          "At the same time, modern pitchers are throwing harder and training more for 
      power and maximum effort. This raises an important question: are pitchers today 
      peaking earlier, declining faster, or experiencing more performance volatility 
      than pitchers from earlier eras?"
        )
      ),
      
      div(
        class = "card",
        h3("Method"),
        p(
          "The analysis uses two main data sources. The modern dataset includes 
      Statcast-era pitcher seasons beginning in 2015. The historical comparison 
      dataset comes from the Lahman Baseball Database and represents a pre-Statcast 
      sample."
        ),
        tags$ul(
          tags$li("Pitcher seasons were filtered to remove very small samples."),
          tags$li("Peak season was defined as each pitcher’s lowest ERA season."),
          tags$li("Mixed-effects models were used to estimate aging curves."),
          tags$li("Natural splines were used so the age-performance relationship could be nonlinear."),
          tags$li("Models controlled for workload, velocity, spin rate, whiff rate, strikeout rate, walk rate, and innings pitched."),
          tags$li("A volatility model measured year-to-year ERA changes."),
          tags$li("A combined model compared Statcast-era and pre-Statcast aging curves.")
        )
      ),
      
      div(
        class = "card",
        h3("Final Results and Conclusion"),
        tags$ul(
          tags$li("Most raw peak seasons occur between ages 25 and 30, with the median peak age around 28."),
          tags$li("The Statcast-era aging curve is fairly stable through the prime years before gradually rising later in a pitcher’s career."),
          tags$li("Higher-velocity pitchers generally show lower predicted ERAs, although velocity does not strongly reshape the entire aging curve by itself."),
          tags$li("Year-to-year ERA volatility changes with age, suggesting that performance stability is an important part of pitcher aging."),
          tags$li("The historical comparison suggests that pre-Statcast pitchers peaked later, while Statcast-era pitchers may peak earlier and decline somewhat faster after age 30."),
          tags$li("Overall, modern pitching intensity may affect not only when pitchers peak, but also how stable and sustainable their performance is over time.")
        )
      ),
      
      div(
        class = "card",
        h3("Why This Matters"),
        p(
          "These findings can help MLB teams think more carefully about pitcher development, 
      workload management, contract timing, and long-term roster planning. If modern 
      pitchers are peaking earlier or becoming less stable later in their careers, teams 
      may need to adjust how they evaluate young pitchers, manage workloads, and project 
      future performance."
        )
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================

server <- function(input, output) {
  
  # ==========================================================
  # DASHBOARD GRAPH OUTPUT
  # ==========================================================
  
  output$mainPlot <- renderPlot({
    
    if (input$plot_type == "Peak Age Distribution") {
      ggplot(pitcher_peaks, aes(x = peak_age)) +
        geom_histogram(binwidth = 1, fill = "#1f5f99", color = "white") +
        labs(
          title = "Distribution of Pitcher Peak Ages",
          subtitle = "Peak defined as each pitcher's lowest ERA season",
          x = "Peak Age",
          y = "Number of Pitchers"
        ) +
        theme_minimal(base_size = 14)
    }
    
    else if (input$plot_type == "Statcast Aging Curve") {
      ggplot(baseline_grid, aes(x = age, y = pred_era)) +
        geom_line(linewidth = 1.4, color = "#1f5f99") +
        labs(
          title = "Estimated Statcast-Era Pitcher Aging Curve",
          subtitle = "Predicted ERA holding workload and pitch characteristics constant",
          x = "Age",
          y = "Predicted ERA"
        ) +
        theme_minimal(base_size = 14)
    }
    
    else if (input$plot_type == "Velocity Aging Curves") {
      ggplot(velocity_grid, aes(x = age, y = pred_era, color = factor(v_z))) +
        geom_line(linewidth = 1.4) +
        labs(
          title = "Predicted Aging Curves at Different Velocity Levels",
          subtitle = "Velocity levels: -1 SD, average, +1 SD",
          x = "Age",
          y = "Predicted ERA",
          color = "Velocity (z)"
        ) +
        theme_minimal(base_size = 14)
    }
    
    else if (input$plot_type == "Volatility Curve") {
      ggplot(vol_grid, aes(x = age_mid, y = pred_vol)) +
        geom_line(linewidth = 1.4, color = "#8b1e3f") +
        labs(
          title = "Predicted Year-to-Year ERA Volatility by Age",
          subtitle = "Volatility measured as absolute ERA change",
          x = "Age",
          y = "Predicted Absolute ERA Change"
        ) +
        theme_minimal(base_size = 14)
    }
    
    else if (input$plot_type == "Pre-Statcast Aging Curve") {
      ggplot(hist_grid, aes(x = age, y = pred_era)) +
        geom_line(linewidth = 1.4, color = "#26734d") +
        labs(
          title = "Pre-Statcast Pitcher Aging Curve",
          subtitle = "Historical comparison sample from Lahman",
          x = "Age",
          y = "Predicted ERA"
        ) +
        theme_minimal(base_size = 14)
    }
    
    else if (input$plot_type == "Era Comparison") {
      ggplot(combined_grid, aes(x = age, y = pred_era, color = era_type)) +
        geom_line(linewidth = 1.4) +
        labs(
          title = "Estimated Pitcher Aging Curves: Pre-Statcast vs Statcast",
          subtitle = "Comparison using common variables across eras",
          x = "Age",
          y = "Predicted ERA",
          color = "Era"
        ) +
        theme_minimal(base_size = 14)
    }
  })
  
  
  # ==========================================================
  # DASHBOARD INTERPRETATION TEXT
  # ==========================================================
  
  output$summaryText <- renderUI({
    
    if (input$plot_type == "Peak Age Distribution") {
      HTML("<p>This histogram shows each pitcher’s raw peak season based on lowest ERA. Most peak seasons cluster in the mid-to-late twenties, supporting the idea that pitcher peak performance generally occurs before the mid-thirties.</p>")
    }
    
    else if (input$plot_type == "Statcast Aging Curve") {
      HTML("<p>This model-based curve estimates predicted ERA across age while holding workload and pitch characteristics constant. The curve is fairly stable through the prime years and rises later in the career, suggesting gradual decline rather than an abrupt collapse.</p>")
    }
    
    else if (input$plot_type == "Velocity Aging Curves") {
      HTML("<p>This graph compares low-, average-, and high-velocity pitchers. Higher-velocity pitchers generally maintain lower predicted ERA values, suggesting that velocity improves performance level even if it does not fully reshape the aging curve.</p>")
    }
    
    else if (input$plot_type == "Volatility Curve") {
      HTML("<p>This curve shows predicted year-to-year ERA instability. Higher values indicate larger changes in ERA from one season to the next. The model suggests that pitcher performance becomes somewhat less stable with age.</p>")
    }
    
    else if (input$plot_type == "Pre-Statcast Aging Curve") {
      HTML("<p>This historical curve provides the pre-Statcast benchmark. It helps compare modern pitcher aging patterns against an earlier period when workload and pitcher usage differed from the current game.</p>")
    }
    
    else if (input$plot_type == "Era Comparison") {
      HTML("<p>This graph directly compares pre-Statcast and Statcast-era aging curves using common variables. The model suggests modern pitchers may peak earlier and decline somewhat faster after age 30.</p>")
    }
  })
  
  
  # ==========================================================
  # PLAYER EXPLORER REACTIVE DATA
  # ==========================================================
  
  selected_player_data <- reactive({
    statcast_clean %>%
      mutate(pitcher_id = as.character(pitcher_id)) %>%
      filter(pitcher_id == input$selected_player) %>%
      arrange(year)
  })
  
  selected_player_peak <- reactive({
    pitcher_peaks %>%
      mutate(pitcher_id = as.character(pitcher_id)) %>%
      filter(pitcher_id == input$selected_player)
  })
  
  
  # ==========================================================
  # PLAYER PROFILE CARD
  # ==========================================================
  
  output$playerProfile <- renderUI({
    pdata <- selected_player_data()
    
    if (nrow(pdata) == 0) {
      return(HTML("<p>No player found.</p>"))
    }
    
    player_name <- pdata$pitcher_name[1]
    player_id <- as.character(pdata$pitcher_id[1])
    img_url <- get_headshot_url(player_id)
    
    div(
      class = "player-card",
      tags$img(
        src = img_url,
        class = "player-img",
        onerror = "this.style.display='none';"
      ),
      div(
        div(class = "player-name", player_name),
        div(class = "player-sub", paste("MLBAM Player ID:", player_id)),
        div(class = "player-sub", paste("Seasons in Dataset:", min(pdata$year), "-", max(pdata$year)))
      )
    )
  })
  
  
  # ==========================================================
  # PLAYER METRIC CARDS
  # ==========================================================
  
  output$playerYearsCard <- renderUI({
    pdata <- selected_player_data()
    
    div(
      class = "metric-card",
      div(class = "metric-value", length(unique(pdata$year))),
      div(class = "metric-label", "Seasons")
    )
  })
  
  output$playerPeakAgeCard <- renderUI({
    ppeak <- selected_player_peak()
    
    div(
      class = "metric-card",
      div(class = "metric-value", ifelse(nrow(ppeak) > 0, ppeak$peak_age[1], "NA")),
      div(class = "metric-label", "Peak Age")
    )
  })
  
  output$playerPeakEraCard <- renderUI({
    ppeak <- selected_player_peak()
    
    div(
      class = "metric-card",
      div(class = "metric-value", ifelse(nrow(ppeak) > 0, round(ppeak$peak_era[1], 2), "NA")),
      div(class = "metric-label", "Best ERA")
    )
  })
  
  output$playerVelocityCard <- renderUI({
    pdata <- selected_player_data()
    
    div(
      class = "metric-card",
      div(class = "metric-value", round(mean(pdata$fastball_avg_speed, na.rm = TRUE), 1)),
      div(class = "metric-label", "Avg Fastball Velo")
    )
  })
  
  
  # ==========================================================
  # PLAYER SEASON PLOT
  # ==========================================================
  
  output$playerSeasonPlot <- renderPlot({
    pdata <- selected_player_data()
    
    ggplot(pdata, aes(x = year, y = p_era)) +
      geom_line(linewidth = 1.2, color = "#1f5f99") +
      geom_point(size = 3, color = "#12355b") +
      labs(
        title = paste("Season-by-Season ERA:", pdata$pitcher_name[1]),
        subtitle = "Lower ERA indicates better run prevention",
        x = "Season",
        y = "ERA"
      ) +
      theme_minimal(base_size = 14)
  })
  
  
  # ==========================================================
  # PLAYER SEASON TABLE
  # ==========================================================
  
  output$playerSeasonTable <- renderTable({
    selected_player_data() %>%
      transmute(
        Year = year,
        Age = age,
        ERA = round(p_era, 2),
        IP = round(p_formatted_ip, 1),
        `K%` = round(k_percent, 1),
        `BB%` = round(bb_percent, 1),
        `Whiff%` = round(whiff_percent, 1),
        `Fastball Velo` = round(fastball_avg_speed, 1),
        `Pitch Count` = pitch_count
      )
  })
}


# ============================================================
# RUN APP
# ============================================================

shinyApp(ui = ui, server = server)
setwd("C:/Users/Alicia/Downloads/SMGT490")
library(rsconnect)
deployApp()

###BloxPlot
library(ggstatsplot)
library(hrbrthemes)
library(ggplot2)

dnam <- read.csv("brs_dataset_incor_laval_norm.txt", sep = "\t")
head(dnam)
row.names(dnam)<- dnam[,1]
dnam<- dnam[,-1]
head(dnam)
dim(dnam)

dnam$Class <- factor(dnam$Class, levels=c("Control", "Possible", "BrS"))

plt <- ggbetweenstats(
  data = dnam,
  x = Class,
  y =Actin, type = "np"
)
plt

plt <- plt + 
  # Add labels and title
  labs(
    x = "Shanghai Score",
    y = "O.D. (Actin)",
    title = "Distribution of Actin O.D. value across Shanghai Score"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "TT Courier New", size = 14, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 18,
      face = "bold",
      color = "black",
      hjust = 0.5
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 14, 
      face = "bold",
      color="black"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

plt

plt <- ggbetweenstats(
  data = dnam,
  x = Class,
  y =Keratin, type = "np"
)
plt

plt <- plt + 
  # Add labels and title
  labs(
    x = "Shanghai Score",
    y = "O.D. (Keratin)",
    title = "Distribution of Keratin O.D. value across Shanghai Score"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "TT Courier New", size = 14, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 18,
      face = "bold",
      color = "black",
      hjust = 0.5
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 14, 
      face = "bold",
      color="black"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

plt
plt <- ggbetweenstats(
  data = dnam,
  x = Class,
  y =Connexin, type = "np"
)
plt

plt <- plt + 
  # Add labels and title
  labs(
    x = "Shanghai Score",
    y = "O.D. (Connexin)",
    title = "Distribution of Connexin O.D. value across Shanghai Score"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "TT Courier New", size = 14, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 18,
      face = "bold",
      color = "black",
      hjust = 0.5
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 14, 
      face = "bold",
      color="black"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

plt
plt <- ggbetweenstats(
  data = dnam,
  x = Class,
  y =Means, type = "np"
)
plt

plt <- plt + 
  # Add labels and title
  labs(
    x = "Shanghai Score",
    y = "O.D. (Means)",
    title = "Distribution of Means O.D. value across Shanghai Score"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "TT Courier New", size = 14, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 18,
      face = "bold",
      color = "black",
      hjust = 0.5
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 14, 
      face = "bold",
      color="black"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

plt

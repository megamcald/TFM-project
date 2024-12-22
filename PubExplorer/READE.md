# PubExplorer Text Mining Analysis
![PubExplorer](PubExplorer/www/logos/logo_appv2.png)

**PubExplorer** es una aplicación web desarrollada en R Shiny para el estudio y análisis de la **enfermedad de Crohn** mediante minería de textos de PubMed. La aplicación permite a los usuarios buscar artículos relevantes, analizar la frecuencia de palabras y genes, visualizar datos en tablas, gráficos y nubes de palabras, y explorar la co-ocurrencia de términos en los abstracts de las publicaciones.

## 📑 Tabla de Contenidos

1. [Descripción](#descripción)
2. [Objetivos](#objetivos)
3. [Instalación](#instalación)
4. [Uso](#uso)
5. [Estructura del Proyecto](#estructura-del-proyecto)
6. [Paquetes y librerías](#paquetes-y-dependencias)
7. [Créditos](#créditos)
9. [Capturas](#capturas-de-pantalla)

## 📝 Descripción

PubExplorer facilita la búsqueda y análisis de artículos relacionados con la enfermedad de Crohn en PubMed. Utilizando técnicas de minería de textos, la aplicación extrae información relevante de los abstracts, analiza la frecuencia de palabras clave y genes asociados, y proporciona visualizaciones para una mejor comprensión de los datos.

## 🎯 Objetivos

El objetivo principal es el desarrollo de una aplicación web interactiva, que use el sistema de búsqueda de PubMed, y permita a los usuarios extraer la información contenida en los abstracts de las publicaciones, aplicando técnicas de minería de textos.

- Facilitar la búsqueda de artículos en PubMed sobre la enfermedad de Crohn.
- Proporcionar estadísticas sobre la frecuencia de palabras clave y genes asociados a la enfermedad.
- Analizar la proximidad de co-ocurrencia entre diferentes términos en los abstracts.
- Generar visualizaciones: nubes de palabras, gráficos de barras y tablas para facilitar el análisis e interpretación de datos.
- Acceso a enlaces directos: PubMed, DOI, QuickGO, UniProt.

## 🛠️ Instalación

### 🎯 Requisitos Previos

- **R** (versión 4.0 o superior)
- **RStudio** (opcional pero recomendado)

### 📦 Paquetes Necesarios

Es necesario tener instalados los siguientes paquetes de R. Se pueden ejecutar con el siguiente código de R:

```r
install.packages(c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "fs",
  "DT",
  "wordcloud",
  "ggplot2",
  "pubmed.mineR",
  "easyPubMed",
  "lsa",
  "tokenizers",
  "viridis"
))
```
### 🔧 Configuración del Proyecto

[Repositorio github](https://github.com/megamcald/TFM-project)🔗

Para clonar el repositorio:

git clone https://github.com/megamcald/TFM-project.git

## 🗂️ Estructura del Proyecto

```r

TFM-project/
├── app.R                     # Archivo principal de la aplicación Shiny
├── README.md                 # README
├── www/
│   ├── styles.css            # Estilos CSS
│   ├── images/               # Capruras aplicación de Shiny
│   ├── logos/				        # Logos
│   │   ├── lupa.png
│   │   ├── logo_appv2.png
│   │   └── logo_uoc.png
│   └── spinnerv2.gif         # GIF del spinner de carga
├── data/
│   ├── raw/
│   │   └── crohns_disease/   # Datos crudos descargados de PubMed
│   └── processed/            # Datos procesados
├── results/                  # Resultados
│   ├── word_frequency/
│   │   ├── word_tokens.txt
│   │   └── word_tokens_barplot.png
│   └── gene_analysis/
│       ├── genes_df.txt
│       └── genes_barplot.png

```
## 📚 Paquetes y librerías

- *Shiny*
	- shiny: Creación de aplicaciones web interactivas.
	- shinydashboard: Creación de paneles de control atractivos.
	- shinyjs: Facilita el uso de JavaScript en Shiny.
- *Manejo de rutas*
	- fs: Manipulación de rutas de archivos y directorios.
- *Visualización*
	DT: Creación de tablas interactivas.
	wordcloud: Generación de nubes de palabras.
	ggplot2: Creación de gráficos estáticos y dinámicos.
- *Acceso PubMed*
	- pubmed.mineR: Acceso y análisis de datos de PubMed.
	- easyPubMed: Consultas y descargas desde PubMed.
- *Minería de textos*
	- lsa: Análisis semántico latente.
	- tokenizers: Tokenización de textos.
- *Otros*
	- viridis: Paletas de colores para visualizaciones.

### 👥 Créditos
Este proyecto fue desarrollado como parte del Trabajo Final del Máster en Bioinformática y Bioestadística en la UOC-UB (Universitat Oberta de Catalunya - Universitat de Barcelona).

Desarrollado por: [Megam Calderón](https://www.linkedin.com/in/megam-calder%C3%B3n/")🔗

### 🖼️ Capturas

## Pantalla principal
![Captura pantalla de inicio PubExplorer](PubExplorer/www/images/vista_inicio.png)


## Pantallas Palabras
![Captura pantalla Palabras-tabla](PubExplorer/www/images/vista_palabras.png)
![Captura pantalla Palabras-gráfico](PubExplorer/www/images/vista_palabras_plot.png)
![Captura pantalla Palabras-nube](PubExplorer/www/images/vista_palabras_cloud.png)

## Pantallas Genes
![Captura pantalla Genes-tabla](PubExplorer/www/images/vista_genes.png)
![Captura pantalla Genes-gráfico](PubExplorer/www/images/vista_genes_plot.png)
![Captura pantalla Genes-nube](PubExplorer/www/images/vista_genes_cloud.png)

## Pantallas Navegación
![Captura pantalla Navegación términos: Gene_symbol](PubExplorer/www/images/vista_navgenes.png)
![Captura pantalla Navegación términos: words](PubExplorer/www/images/vista_navpalabras.png)

## Pantalla Co-ocurrencia términos
![Captura pantalla Búsqueda co-ocurrencia términos](PubExplorer/www/images/vista_coocurrencia.png)

## Pantalla Acerca de
![Captura pantalla acerca de](PubExplorer/www/images/vista_acercade.png)

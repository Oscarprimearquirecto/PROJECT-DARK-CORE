
# PROJECT DARK-CORE V2.0

> **Pipeline bioinformático in silico de bajo cómputo para la caracterización tensorial, modelado 3D y estimación de afinidad en microproteínas (ncORFs) del Proteoma Oscuro.**

---

## 📌 Descripción General

**PROJECT DARK-CORE** es una arquitectura analítica ligera (*headless*) diseñada para identificar, priorizar y simular microproteínas codificadas por marcos de lectura abiertos no canónicos (ncORFs). 

A diferencia de las plataformas tradicionales de biología estructural que exigen infraestructura HPC masiva, este motor utiliza mapeo tensorial de aminoácidos, cálculo de entropía de gradientes y modelado geométrico aproximado para ejecutar cribados (*screening*) de alta velocidad en entornos de cómputo restringidos (incluyendo arquitecturas ARM64 y entornos móviles).

---

## 📐 Fundamento Matemático y Módulos

El pipeline procesa la secuencia primaria de aminoácidos a través de 4 módulos integrados:

1. **Vectorización Tensorial & Gradientes:** Mapeo de la secuencia a la escala de hidrofobicidad Kyte-Doolittle $V \in \mathbb{R}^N$ y cálculo del gradiente discreto $\nabla V$.
2. **Structural Viability Score (SVR) & Entropía:** Integración de la hidrofobicidad media, la desviación estándar del gradiente y la entropía de Shannon ($H$) sobre la distribución del gradiente:
   $$SVR = 0.4 \cdot \bar{V} + 0.3 \cdot \sigma(\nabla V) + 0.3 \cdot H(\nabla V)$$
3. **Generación Geométrica Backbone (.PDB):** Construcción aproximada de coordenadas tridimensionales en formato estándar PDB ATOM para evaluación de estructura secundaria.
4. **Docking Estimado ($\Delta G$):** Cálculo de energía libre de unión aproximada en kcal/mol frente a dianas estructurales.

---

## 💻 Uso Rápido

```bash
# Ejecutar el orquestador principal
python dark_core_v2.py

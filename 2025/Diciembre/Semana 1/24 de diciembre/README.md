# Quadratic Equation Analyzer

## Project Description

This project implements an algorithm that receives a quadratic equation of the form:

ax^2 + bx + c = 0

and performs a mathematical analysis to determine its fundamental properties.

---

## Features

The algorithm computes and returns:

- Discriminant of the equation
- Type of roots:
  - Two distinct real roots
  - One real (double) root
  - Complex roots
- Real roots, if they exist in ℝ
- Factorization, only when possible over the real numbers

---

## Mathematical Criteria

The discriminant is calculated as:

Δ = b² - 4ac

Root classification:

- Δ > 0 → Two distinct real roots  
- Δ = 0 → One real double root  
- Δ < 0 → Complex roots (no factorization in ℝ)

---

## Constraints

- Factorization is performed only when the equation has real roots.
- If the roots are complex, the algorithm reports the result but does not attempt factorization.

---

## Project Goals

- Strengthen the relationship between algebra and programming
- Practice clean, modular, and testable code
- Build a reusable mathematical analysis tool
- Apply proper validation and error handling

---

## Status

In development — academic and personal use

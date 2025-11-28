---
# An instance of the Experience widget.
# Documentation: https://wowchemy.com/docs/page-builder/
widget: experience

# This file represents a page section.
headless: true

# Order that this section appears on the page.
weight: 40

title: Experience
subtitle:

# Date format for experience
#   Refer to https://wowchemy.com/docs/customization/#date-format
date_format: Jan 2006

# Experiences.
#   Add/remove as many `experience` items below as you like.
#   Required fields are `title`, `company`, and `date_start`.
#   Leave `date_end` empty if it's your current employer.
#   Begin multi-line descriptions with YAML's `|2-` multi-line prefix.
experience:
  - title: Observational Research Manager
    company: Amgen
    company_url: 'https://www.amgen.com'
    company_logo: Amgen
    location: Thousand Oaks, CA (remote in Denver, CO)
    date_start: '2024-07-15'
    date_end: ''
    description: |2-
        Major responsibilities include:
        
        * Regulatory Real-World Evidence (RWE) & PMR Leadership
          * Led design and authorship of FDA-required studies, including a post-marketing requirement final report, NCO study protocols/SAPs, and algorithms validating 3-point MACE and CV mortality using Medicare–REGARDS linked data
          * Developed and led the clean-room gating framework with quantitative assessments of potential unmeasured confounding to support upcoming comparative safety work
          * Conducted multiple ad-hoc analyses integrating claims and biometric EHR data through multiple imputation and transportability methods, demonstrating negligible unmeasured confounding due to biometric covariates in upcoming claims-only safety study
          * Ran simulation studies evaluating bias from nondifferential outcome misclassification and applied validation-based adjustments to recover true effects
          * Authored several regulatory briefing documents and RTQ responses, ensuring methodological clarity and alignment with FDA expectations
        
        * Additional Regulatory, Methodological, and Cross-Functional Contributions
          * Co-led protocol development for the Effectiveness NCO study, aligning core design elements with the Safety NCO study and adding two new negative control outcomes
          * Supported global regulatory needs through regulator-mandated safety analyses across U.S. and Japan datasets, identifying cross-regional methodological efficiencies
          * Contributed to Japan PMO/MOP work through study feasibility assessments and creation of an R Shiny app for automated power/sample size estimation across multiple testing frameworks
          * Authored key sections of two briefing documents for Japanese biopharma regulator, clarifying misclassification bias and differential attrition concerns in regulator-conducted studies
          * Coauthored two manuscripts on postoperative outcomes among women with osteoporosis and supported broader publication strategy with bone and publication teams
          * Evaluated PCORNet network performance to optimize site selection for future claims–EHR linked analyses
          * Served as a reviewer for Amgen’s internal protocol governance body

        
  - title: Senior epidemiologist
    company: Target RWE
    company_url: 'https://www.targetrwe.com'
    company_logo: targetrwe
    location: Durham, NC (remote in NY)
    date_start: '2022-07-01'
    date_end: '2024-07-04'
    description: |2-
        Responsibilities included:
        
        * Designing pharmacoepidemiological hepatology, dermatology, and oncology studies with claims and electronic health record data, including studies that integrate both types of data
        * Implementing advanced causal inference methods to surmount common problems in pharmacoepidemiological analyses
          * Clone censor weight approach to address immortal time bias due to the start of follow-up and treatment initiation not coinciding
          * Inverse probability of censoring weight approach for informative right censoring
          * Prevalent new user design for non-contemporaneous marketing in comparative effectiveness studies
        * Writing proposals and developing presentations to help Target RWE win new contracts and continue existing ones

  - title: Doctoral candidate
    company: Columbia University
    company_url: 'https://www.publichealth.columbia.edu/'
    company_logo: mailman
    location: New York, NY
    date_start: '2017-08-01'
    date_end: '2023-03-28'
    description: |2-
        * Completed doctoral epidemiology training and dissertation on nicotine vaping's unintended consequences
        * Performed teaching assistant duties for numerous epidemiology and statistics courses

design:
  columns: '2'
---

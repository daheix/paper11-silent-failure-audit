# Cover Letter Draft (投稿前人工定稿; v2 — 反映 M1 源码级重定性 + COMPEL 定位)

Dear Editor,

We submit "Silent Failures in Frequency-Domain Finite-Element Analysis of
Cage Induction Motors: A Taxonomy, Forensic Case Studies, and an Audit
Protocol" for consideration in COMPEL — The International Journal for
Computation and Mathematics in Electrical and Electronic Engineering.

Using the open-source GetDP 4.0.0 toolchain (PETSc 3.25.5; MUMPS 5.8.2;
serial SuperLU), we study a class of simulation failures in which the
solver returns, no warning is raised, and yet the reported physical
quantity is wrong. We contribute:

(1) A taxonomy of three silent-failure modes in coupled cage-motor FE
    models, each documented as a forensic case with archived evidence:
    (A) a MUMPS workspace-budget abort (INFO(1) = -19, ICNTL(23),
    identified at source level in MUMPS 5.8.2 dfac_b.F:283) whose all-NaN
    field is silently propagated downstream on a system that an
    independent serial-SuperLU factorization solves in 26 s (residual
    2.28e-5); (B) a "half-field" solution at an intermediate geometric
    scale, physically wrong by 42 percent with no numerical diagnostic;
    (C) a 90-degree sector reduction that violates a uniform-sign
    feasibility condition on harmonic commensurability, inflating the
    air-gap fundamental by 1.44x and torque by 1.54x while KCL remains
    at display precision zero.

(2) A formalization of the uniform-sign feasibility condition for sector
    reduction (previously folklore), with an executable check
    (asm_spec.sector_feasibility) that delineates admissible sector
    angles — 90 degrees is inadmissible for the 12-slot/4-pole family
    while 180 degrees passes.

(3) A five-step audit protocol (completeness, numerical hygiene,
    cross-solver verdict, scale-band regression, uniform-sign pre-check)
    with an executable reference implementation whose doctests exercise
    both the archived forensic data and synthetic fixtures for the
    NaN-propagation and cross-solver-mismatch paths.

All quantitative claims trace to archived raw data bundled with the
manuscript sources; a data map (data/raw/README.md) links every figure
and number to its artefact. We believe the topic fits COMPEL's scope in
computational electromagnetics and numerical methods of electrical
machines.

The manuscript is not under consideration elsewhere. All authors have
approved the submission. Data and code availability: repository archive
bundled with the submission.

Suggested reviewers: (to be completed by the author at submission time)

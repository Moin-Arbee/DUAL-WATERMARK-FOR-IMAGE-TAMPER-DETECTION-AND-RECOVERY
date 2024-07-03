# DUAL-WATERMARK-FOR-IMAGE-TAMPER-DETECTION-AND-RECOVERY
Dual Watermark for Image Tamper Detection and Recovery
This repository contains the implementation of a dual watermarking scheme designed for robust image tamper detection and recovery. The project, developed in MATLAB, embeds a 12-bit watermark into the least significant bits (LSBs) of each pixel in partner blocks to ensure the integrity and authenticity of digital images.

Features
Watermark Embedding:

A 12-bit watermark is embedded into the three LSBs of each pixel in the corresponding partner blocks.
The same watermark is embedded into two blocks to enhance robustness and redundancy.
Tamper Detection:

A four-level detection mechanism is employed to identify tampered regions.
The multi-level approach aims to reduce false alarms and accurately localize tampered areas.
Tamper Recovery:

A two-stage recovery process is implemented to restore tampered blocks.
The recovery process ensures the image's integrity and authenticity post-tampering.
Modules
Watermark Embedding:

Embeds the computed 12-bit watermark into the LSBs of the partner blocks.
Generates the watermarked image from the host image.
Tamper Detection with Localization:

Detects tampered regions through a four-level detection process.
Reduces false positives and accurately localizes tampered areas.
Recovery of Tampered Blocks:

Restores the tampered blocks using a two-stage recovery mechanism.
Ensures the integrity and authenticity of the recovered image.

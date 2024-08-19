# NullControllability_MBSA_DroneReliability
This repository contains the essential codes and system architectures used in the research paper titled "Integrating Null Controllability and Model-Based Safety Assessment for Enhanced Reliability in Drone Design." The repository is divided into two main sections:

1. MATLAB Repository
This section includes two key MATLAB R2022b files for performing an ACAI-based controllability assessment of a selected UAV design:
- InCtrb_AltaRica_rev2.m: This is the main code that provides a GUI for selecting the UAV design, AltaRica model, and controllability hypothesis. It integrates the controllability results with the chosen AltaRica 3.0 model to generate an updated model.
- DB_designs.m: This file contains the database of UAV designs, which is accessed by InCtrb_AltaRica_rev2.m to perform the controllability assessment and integration.

2. System Analyst Repository
This section contains the UAV system architecture models created in System Analyst 1.3, corresponding to the designs from iterations 1 to 6 for Scenario 1 and iterations 1 to 5 for Scenario 2, as discussed in the paper. Each model illustrates the progression of design changes across the iterations.

# Step-by-Step Guide
Step 1 - Create or Open a System Model:
Use System Analyst 1.3 to create a new system model or open an existing one from the System Analyst repository.

Step 2 - Compile the System Model:
In System Analyst 1.3, compile the system model into an AltaRica 3.0 model, generating a .alt file.

Step 3 - Run the MATLAB Script:
Open MATLAB R2022b and run InCtrb_AltaRica_rev2.m.

Step 4 - Use the GUI to Configure the Assessment:
Via the GUI, select your UAV design. Choose the maximum number of simultaneous failures to be considered. Decide whether to include the Yaw axis controllability requirement. Select the AltaRica 3.0 model you created in Step 2.

Step 5 - Generate the New AltaRica 3.0 Model:
Choose a name for the new AltaRica 3.0 model that will integrate the selected design with the controllability assessment.

Step 6 - Compile and Compute Reliability Metrics:
Use AltaRicaWizard 1.2.0 to compile the new AltaRica 3.0 model created in Step 5 and compute the reliability metrics.

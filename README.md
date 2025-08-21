<div align="center">

# Tracking of discrete-time unmodeled reference signals in robotic manipulators using output regulation theory and high-gain observers  
### Manuscript ID: 9931    
</div>
Ricardo Tapia-Herrera¹, Tonatiuh Hernández-Cortés², Alexis Rojas-Ruíz³, Beatriz A. Jaime-Fonseca⁴, Jesús A. Meda-Campaña³  

¹ SECIHTI-IPN, Ciudad de México, México  
² Universidad Politécnica de Pachuca, Zempoala, México  
³ Instituto Politécnico Nacional, SEPI ESIME Zacatenco, Ciudad de México, México  
⁴ Instituto Politécnico Nacional, ESIME Zacatenco, Ciudad de México, México  

---

## 🛠 Requirements

To run the programs you need:

- **MATLAB/Simulink R2024a** or later  
- **Multibody Simulation Toolbox** (required to load the robot’s multibody model)
- **Robotics Toolbox** (required to generate the linear interpolation in Case II)  

> ⚠️ The Simulink files are linked to MATLAB functions that contain the parameters of the controller and plotting functions. We recommend downloading the **complete directory** to ensure all dependencies are included.

---

## Included files
📂 Main simulation files 
| Simulink File           | Related Figure(s)         | Description                                                                 |
|--------------------------|---------------------------|-----------------------------------------------------------------------------|
| `HGO_robot_CASE_I.slx`   | Fig. 4, Fig. 5           | Executes the simulation for Case I: Tracking of a chaotic dynamical system  |
| `HGO_robot_CASE_II.slx`  | Fig. 7, Fig. 8, Fig. 9, Fig. 10 | Executes the simulation for Case II: Tracking of a prescribed workspace trajectory |
| `HGO_robot_CASE_III.slx` | Fig. 11, Fig. 12, Fig. 13, Fig. 14 | Executes the simulation for Case III: Tracking of a reference signal generated in real time |  


📂 Additional scripts

| MATLAB Script                  | Description                                                                 |
|--------------------------------|-----------------------------------------------------------------------------|
| `HGO_parameters_Case_I.m`      | Computes the reference signal for the Lorenz attractor in the robot's workspace. Also calculates the inverse kinematics, stabilizer gains, HGO parameters, and solves Isidori's equations. |
| `HGO_parameters_Case_II.m`     | Computes the reference signal for a trajectory with five linear segments. Also calculates the inverse kinematics, stabilizer gains, HGO parameters, and solves Isidori's equations. |
| `HGO_parameters_Case_III.m`    | Computes the reference signal for real-time tracking. Also calculates inverse kinematics, stabilizer gains, HGO parameters, and solves Isidori's equations. |
| `Lorenz_attractor_simulink.slx`| Required by `HGO_parameters_Case_I.m` to compute the states of the Lorenz attractor for reference tracking. |
| `Lorenz_attractor_fnc.m`       | Contains the model of the Lorenz attractor and is used by the Simulink file `Lorenz_attractor_simulink.slx`. |
| `results.m`                    | Generates the plots presented in the article. This script runs automatically once a main Simulink file completes the simulation. |  
| `px_final.mat` and 'py_final'  | These files contain the adquired values in the robot's workspace (Px, Py) for Case III. They are automatically loaded when 'HGO_robot_CASE_III.slx is executed'|  




## 📊 Results

- Once the Simulink file finishes execution, the **plots of the results** will appear automatically.  
- The figures correspond to the tracking performance and output regulation behavior discussed in the paper.  

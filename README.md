# fNIRS–Stimulus Correlation Analysis

This project demonstrates that **the fNIRS signal is correlated with the stimulus** using a public dataset and standard tools (MNE-Python). It includes a full preprocessing pipeline, epoch-based evoked responses, and explicit correlation measures between the haemodynamic signal and task condition.

---

## Dataset

We use the **MNE fNIRS motor dataset** ([MNE documentation](https://mne.tools/stable/documentation/datasets.html#fnirs-motor)), a well-documented public dataset recorded at Macquarie University.

- **Design**: Block design; 5 s blocks per trial.
- **Conditions**:
  - **Tapping/Left**: participant taps thumb to fingers on the left hand.
  - **Tapping/Right**: same with the right hand.
  - **Control**: no motor task.
- **Trials**: 30 per condition (90 total); trigger codes 1 = Control, 2 = Tapping/Left, 3 = Tapping/Right (code 15 = experiment start/end, excluded).
- **Montage**: Optodes over motor cortex; 56 fNIRS (CW amplitude) channels, ~7.8 Hz sampling.

The data are downloaded automatically on first run via `mne.datasets.fnirs_motor.data_path()`.

---

## Analysis Pipeline

### 1. Preprocessing (standard fNIRS in MNE)

1. **Annotations**: Trigger codes mapped to condition names; stimulus duration set to 5 s; experiment start/end (15) removed.
2. **Channel selection**: Short source–detector pairs (< 1 cm) removed (insensitive to cortical signal).
3. **Intensity → optical density**: `mne.preprocessing.nirs.optical_density()`.
4. **Optical density → haemoglobin**: Modified Beer–Lambert law (`beer_lambert_law`, ppf=0.1) → relative HbO and HbR.
5. **Bandpass filter**: 0.05–0.7 Hz to retain haemodynamic range and attenuate cardiac (~1 Hz) and drift.

### 2. Epoching

- Epochs from −2 s to +12 s around each stimulus onset.
- Baseline: pre-stimulus (None to 0).
- Rejection: epochs with HbO > 80 µM rejected as artifact.

### 3. Assessing correlation with stimulus

Three complementary approaches:

1. **Evoked response (task vs control)**  
   Average HbO/HbR across trials and channels for “Tapping” (Left + Right) vs “Control”. Compare mean amplitude in the typical haemodynamic peak window (4–10 s).

2. **Continuous stimulus regressor**  
   Binary regressor: 1 during any tapping block, 0 otherwise, aligned to the continuous mean HbO time series. **Pearson correlation** between this regressor and the (mean across channels) HbO signal.

3. **Trial-level correlation**  
   For each epoch, compute mean HbO in the 4–10 s window. **Pearson correlation** between condition (1 = Tapping, 0 = Control) and this peak-HbO value.

---

## How to run

```bash
pip install -r requirements.txt   # or: python3 -m pip install -r requirements.txt
python fnirs_stimulus_correlation_analysis.py
```

On first run, the script will download the dataset (~18 MB) to `~/mne_data`. Figures are saved in `results/`.

---

## Results

### Evoked responses

- **Mean HbO in 4–10 s** (averaged across channels): **Tapping ≈ 6.83 µM**, **Control ≈ −0.87 µM** → difference **~7.7 µM**.
- HbO increases during tapping and is flat/slightly negative during control, consistent with motor-cortex activation.
- HbR shows the expected decrease (not shown in the summary stats but visible in the evoked plot).

### Stimulus–signal correlation

- **Stimulus regressor vs mean HbO (continuous)**: **r ≈ 0.044**, **p ≈ 2.7e−11**.  
  The correlation is small in magnitude (long, noisy time series) but **highly significant**, indicating that the binary “tapping on/off” regressor aligns with HbO fluctuations.

- **Trial-level (condition vs peak HbO in 4–10 s)**: **r ≈ 0.41**, **p ≈ 1.2e−04**.  
  Condition strongly predicts peak HbO: tapping trials show higher HbO in the haemodynamic window than control trials.

### Figures (in `results/`)

- **`evoked_tapping_vs_control.png`**: Evoked HbO and HbR for Tapping vs Control; shaded region = 4–10 s window.
- **`stimulus_signal_timeseries.png`**: First 200 s of stimulus regressor and mean HbO (alignment and correlation).
- **`epoch_condition_vs_peak_hbo.png`**: Scatter of condition (0/1) vs peak HbO per epoch.

---

## Discussion

- **Stimulus–signal relationship**: The analysis shows a clear link between the fNIRS signal and the stimulus:
  - Evoked HbO is larger for Tapping than for Control in the 4–10 s window.
  - The continuous stimulus regressor is significantly correlated with mean HbO.
  - At the trial level, condition predicts peak HbO (moderate r, significant p).

- **Interpretation**: The increase in HbO over motor cortex during finger tapping is consistent with neurovascular coupling: motor activity increases local metabolism and blood flow, and fNIRS detects the resulting haemodynamic response. The control condition shows no such rise.

- **Limitations**: Single subject; one public dataset; no correction for multiple comparisons across channels. The continuous r is small because the haemodynamic response is slow and the regressor is binary; trial-level and evoked analyses are better suited to show the effect.

- **Tools**: All steps use **MNE-Python** (reading, preprocessing, epoching, averaging) and **SciPy** (correlation). No custom or proprietary tools are required.

---

## References and links

- MNE-Python: [mne.tools](https://mne.tools/)
- MNE fNIRS preprocessing tutorial: [Preprocessing fNIRS data](https://mne.tools/stable/auto_tutorials/preprocessing/70_fnirs_processing.html)
- fNIRS public data (NITRC): [fnirsdata](https://www.nitrc.org/projects/fnirsdata/)
- Example BIDS fNIRS (e.g. finger tapping): [BIDS-NIRS-Tapping](https://github.com/rob-luke/BIDS-NIRS-Tapping)

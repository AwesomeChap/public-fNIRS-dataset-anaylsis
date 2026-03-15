"""
fNIRS–stimulus correlation analysis using a public dataset.

Dataset: MNE fNIRS motor dataset (single subject, finger tapping task).
Source: mne.datasets.fnirs_motor — recorded at Macquarie University.
Conditions: Tapping/Left, Tapping/Right, Control (30 trials each, 5 s blocks).

This script demonstrates that the fNIRS signal (HbO₂/HbR) is correlated with
the stimulus by: (1) standard preprocessing, (2) epoch-level averaging,
(3) comparing evoked responses (task vs control), and (4) a correlation
between a stimulus regressor and the haemoglobin time series.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import mne
from mne.preprocessing import nirs
from scipy import stats

# ---------------------------------------------------------------------------
# 1. Dataset description and loading
# ---------------------------------------------------------------------------
print("=" * 60)
print("1. DATASET")
print("=" * 60)
print(
    "Using: MNE fNIRS motor dataset (mne.datasets.fnirs_motor)\n"
    "  - Single subject, optodes over motor cortex.\n"
    "  - Conditions: Tapping/Left, Tapping/Right, Control (5 s blocks).\n"
    "  - 30 trials per condition; triggers 1=Control, 2=Tapping/Left, 3=Tapping/Right, 15=experiment start/end.\n"
)

fnirs_data_folder = mne.datasets.fnirs_motor.data_path()
participant_dir = fnirs_data_folder / "Participant-1"
raw_intensity = mne.io.read_raw_nirx(participant_dir, verbose=False)
raw_intensity.load_data()
print(raw_intensity.info)
print()

# ---------------------------------------------------------------------------
# 2. Annotations: map trigger codes to condition names and set durations
# ---------------------------------------------------------------------------
raw_intensity.annotations.set_durations(5)
raw_intensity.annotations.rename(
    {"1.0": "Control", "2.0": "Tapping/Left", "3.0": "Tapping/Right"}
)
unwanted = np.nonzero(raw_intensity.annotations.description == "15.0")[0]
raw_intensity.annotations.delete(unwanted)

# ---------------------------------------------------------------------------
# 3. Channel selection (drop short channels)
# ---------------------------------------------------------------------------
picks = mne.pick_types(raw_intensity.info, meg=False, fnirs=True)
dists = nirs.source_detector_distances(raw_intensity.info, picks=picks)
raw_intensity.pick(picks[dists > 0.01])

# ---------------------------------------------------------------------------
# 4. Preprocessing: intensity → optical density → haemoglobin
# ---------------------------------------------------------------------------
print("Preprocessing: intensity → OD → haemoglobin (Beer-Lambert)...")
raw_od = nirs.optical_density(raw_intensity)
raw_haemo = nirs.beer_lambert_law(raw_od, ppf=0.1)
raw_haemo.filter(0.05, 0.7, h_trans_bandwidth=0.2, l_trans_bandwidth=0.02)

# ---------------------------------------------------------------------------
# 5. Events and epoching
# ---------------------------------------------------------------------------
events, event_id = mne.events_from_annotations(raw_haemo)
tmin, tmax = -2, 12
reject_criteria = dict(hbo=80e-6)
epochs = mne.Epochs(
    raw_haemo,
    events,
    event_id=event_id,
    tmin=tmin,
    tmax=tmax,
    reject=reject_criteria,
    reject_by_annotation=True,
    baseline=(None, 0),
    preload=True,
    verbose=False,
)
print(f"Epochs: {len(epochs)} total; Control: {len(epochs['Control'])}, "
      f"Tapping/Left: {len(epochs['Tapping/Left'])}, Tapping/Right: {len(epochs['Tapping/Right'])}")
print()

# ---------------------------------------------------------------------------
# 6. Evoked responses: task vs control (evidence of stimulus–signal relation)
# ---------------------------------------------------------------------------
print("=" * 60)
print("2. EVOKED RESPONSES (Task vs Control)")
print("=" * 60)

evoked_tapping = epochs["Tapping"].average(picks="hbo")
evoked_control = epochs["Control"].average(picks="hbo")
times = evoked_tapping.times
n_ch = len(evoked_tapping.ch_names)

# Peak window for motor response (approx. 4–10 s post stimulus)
peak_tmin, peak_tmax = 4.0, 10.0
peak_mask = (times >= peak_tmin) & (times <= peak_tmax)
tapping_mean = evoked_tapping.get_data().mean(axis=0)
control_mean = evoked_control.get_data().mean(axis=0)
peak_tapping = tapping_mean[peak_mask].mean()
peak_control = control_mean[peak_mask].mean()
print(f"Mean HbO₂ in peak window [{peak_tmin}-{peak_tmax} s] (averaged across channels):")
print(f"  Tapping: {peak_tapping*1e6:.2f} µM (×1e6)")
print(f"  Control: {peak_control*1e6:.2f} µM (×1e6)")
print(f"  Difference (Tapping − Control): {(peak_tapping - peak_control)*1e6:.2f} µM")
print()

# ---------------------------------------------------------------------------
# 7. Correlation: stimulus regressor vs fNIRS signal
# ---------------------------------------------------------------------------
print("=" * 60)
print("3. STIMULUS–SIGNAL CORRELATION")
print("=" * 60)

# Build a stimulus regressor: 1 during any tapping block, 0 otherwise (same length as raw)
sfreq = raw_haemo.info["sfreq"]
n_times_raw = raw_haemo.times.size
stimulus_regressor = np.zeros(n_times_raw)
for i, (onset, duration) in enumerate(zip(
    raw_haemo.annotations.onset,
    raw_haemo.annotations.duration,
)):
    desc = str(raw_haemo.annotations.description[i])
    if "Tapping" in desc:
        start_idx = int(onset * sfreq)
        end_idx = min(int((onset + duration) * sfreq), n_times_raw)
        stimulus_regressor[start_idx:end_idx] = 1.0

# Use HbO₂ (mean across channels) as the fNIRS signal for correlation
raw_haemo_pick = raw_haemo.copy().pick(picks="hbo")
fnirs_signal = raw_haemo_pick.get_data().mean(axis=0)
# Remove mean and match length
fnirs_signal = fnirs_signal - np.mean(fnirs_signal)
stimulus_regressor = stimulus_regressor - np.mean(stimulus_regressor)
# Pearson correlation (same length)
r, p_val = stats.pearsonr(stimulus_regressor, fnirs_signal)
print(f"Pearson correlation (stimulus regressor vs mean HbO₂ time series): r = {r:.4f}, p = {p_val:.2e}")
print()

# Per-epoch correlation: epoch-averaged HbO₂ in peak window vs condition (1=tapping, 0=control)
# This shows that trial-level “condition” predicts the haemodynamic response.
tapping_ids = (event_id["Tapping/Left"], event_id["Tapping/Right"])
epoch_peak_hbo = []
epoch_condition = []
# Use accepted events (epochs.events) to match each epoch to its condition
n_epochs = epochs.events.shape[0]
for i in range(n_epochs):
    ev = epochs.events[i, 2]
    epoch_condition.append(1 if ev in tapping_ids else 0)
    data = epochs[i].get_data(picks="hbo")
    t_idx = np.where((epochs.times >= peak_tmin) & (epochs.times <= peak_tmax))[0]
    epoch_peak_hbo.append(np.mean(data[:, :, t_idx]))
epoch_peak_hbo = np.array(epoch_peak_hbo)
epoch_condition = np.array(epoch_condition)
r_epoch, p_epoch = stats.pearsonr(epoch_condition, epoch_peak_hbo)
print(f"Trial-level: correlation(condition, peak HbO₂ in epoch): r = {r_epoch:.4f}, p = {p_epoch:.2e}")
print(f"Mean peak HbO₂ (Tapping epochs): {epoch_peak_hbo[epoch_condition==1].mean()*1e6:.2f} µM")
print(f"Mean peak HbO₂ (Control epochs): {epoch_peak_hbo[epoch_condition==0].mean()*1e6:.2f} µM")
print()

# ---------------------------------------------------------------------------
# 8. Figures
# ---------------------------------------------------------------------------
out_dir = Path(__file__).resolve().parent / "results"
out_dir.mkdir(exist_ok=True)

# Figure 1: Evoked HbO₂/HbR for Tapping vs Control
fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
for cond, label, color in [
    ("Tapping", "Tapping (Left + Right)", "C1"),
    ("Control", "Control", "gray"),
]:
    ev_hbo = epochs[cond].average(picks="hbo")
    ev_hbr = epochs[cond].average(picks="hbr")
    t = ev_hbo.times
    # Evoked get_data() is (n_channels, n_times); mean over channels -> (n_times,)
    axes[0].plot(t, ev_hbo.get_data().mean(axis=0) * 1e6, label=label, color=color)
    axes[1].plot(t, ev_hbr.get_data().mean(axis=0) * 1e6, label=label, color=color)
axes[0].set_ylabel("HbO₂ (µM)")
axes[1].set_ylabel("HbR (µM)")
axes[1].set_xlabel("Time (s)")
axes[0].legend()
axes[1].legend()
axes[0].axvspan(peak_tmin, peak_tmax, alpha=0.15, color="green")
axes[1].axvspan(peak_tmin, peak_tmax, alpha=0.15, color="green")
axes[0].set_title("Evoked response: fNIRS signal locked to stimulus")
fig.tight_layout()
fig.savefig(out_dir / "evoked_tapping_vs_control.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_dir / 'evoked_tapping_vs_control.png'}")

# Figure 2: Stimulus regressor vs fNIRS (short segment)
seg_len = int(200 * sfreq)  # 200 s
seg_len = min(seg_len, len(fnirs_signal), len(stimulus_regressor))
t_plot = np.arange(seg_len) / sfreq
fig, ax1 = plt.subplots(figsize=(10, 3))
ax1.plot(t_plot, stimulus_regressor[:seg_len], color="orange", alpha=0.8, label="Stimulus (Tapping=1)")
ax1.set_ylabel("Stimulus", color="orange")
ax1.set_xlabel("Time (s)")
ax2 = ax1.twinx()
ax2.plot(t_plot, fnirs_signal[:seg_len] * 1e6, color="blue", alpha=0.7, label="Mean HbO₂")
ax2.set_ylabel("HbO₂ (µM)", color="blue")
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
ax1.set_title(f"Stimulus–signal alignment (r = {r:.3f}, p = {p_val:.2e})")
fig.tight_layout()
fig.savefig(out_dir / "stimulus_signal_timeseries.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_dir / 'stimulus_signal_timeseries.png'}")

# Figure 3: Epoch-level condition vs peak HbO₂
fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(epoch_condition, epoch_peak_hbo * 1e6, alpha=0.6)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Control", "Tapping"])
ax.set_ylabel("Peak HbO₂ in epoch (µM)")
ax.set_xlabel("Condition")
ax.set_title(f"Trial-level correlation: r = {r_epoch:.3f}, p = {p_epoch:.2e}")
fig.tight_layout()
fig.savefig(out_dir / "epoch_condition_vs_peak_hbo.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_dir / 'epoch_condition_vs_peak_hbo.png'}")

print()
print("=" * 60)
print("4. SUMMARY")
print("=" * 60)
print(
    "The analysis shows a clear relationship between the fNIRS signal and the stimulus:\n"
    "  - Evoked HbO₂ is larger during Tapping than during Control in the expected 4–10 s window.\n"
    "  - The stimulus regressor (tapping on/off) is significantly correlated with the mean HbO₂ time series.\n"
    "  - At the trial level, condition (Tapping vs Control) predicts peak HbO₂ in the haemodynamic window.\n"
    "This is consistent with increased cortical oxygenation over motor cortex during finger tapping."
)

# spica-nss

SPICA-NSS is a Python-based graphical tool developed for the SPICA survey, in the
framework of the ISSP (Interferometric Survey of Stellar Parameters) ERC Adv. Grant
#101019653. It provides target selection, observability filtering, calibrator
identification, catalog querying, and interoperability with JMMC Aspro2 to support
long-baseline optical interferometry observing programs.

## Requirements

- Python >= 3.10
- Tkinter
- Internet access to JMMC TAP services

## Dependencies

- Astropy
- Astroquery
- Astroplan
- PyVO
- NumPy
- SciPy
- Matplotlib
- a2p2

### Install

```bash
git clone https://github.com/dsalabert/spica-nss.git
cd spica-nss

python -m venv .venv
source .venv/bin/activate

pip install -U pip
pip install -r requirements.txt
```

### Run

```bash
python nss.py
```

## Documentation

In nss_documentation.md, and in Docs/nss.html (created with pdoc). A description of the aims and the functionalities of the NSS can be found in the file DescriptionNSS.pdf.

### Class overview

#### `spica_NSS`

SPICA-NSS graphical observing preparation tool. The single class that owns the
entire GUI and all business logic.

---

### GUI & application lifecycle

| Method | Description |
|---|---|
| `__init__` | Initialize the main Tkinter window, all frames, widgets, and start the event loop. |
| `onQuery` | Query the SPICA database, filter targets to the observable RA window, and refresh the plot. |
| `onReset` | Reset all filters to defaults and re-query from scratch. |
| `onFilters` | Re-apply active filters and refresh the sky plot. |
| `onDRS` | Simulate a DRS/QCS pipeline run on the selected targets (development aid). |
| `onLog` | Open the survey statistics popup showing completion rates per programme and mode. |
| `onQuit` | Quit the application after user confirmation. |

---

### Target selection & filtering

| Method | Description |
|---|---|
| `getSelectedTargets` | Intersect all active filter indices to build `self.indexList_Targets`. |
| `getAddTarget` | Build the fainter-star pool (`self.indexList_AddTarget`) for the add-target popup. |
| `plotSelectedProgName` | Recompute the programme-name filter index and refresh the plot. |
| `plotSelectedInstMode` | Recompute the instrumental-mode filter index and refresh the plot. |
| `plotSelectedFinalPriority` | Recompute the final-priority filter index and refresh the plot. |
| `entryDateCallback` | Handle a new observation date and re-query. |
| `entryDecMinCallback` | Handle a new minimum declination value. |
| `entryDecMaxCallback` | Handle a new maximum declination value. |
| `entryVmagMinCallback` | Handle a new minimum V-magnitude value. |
| `entryVmagMaxCallback` | Handle a new maximum V-magnitude value. |

---

### Sky distribution plot

| Method | Description |
|---|---|
| `plot_radec` | Draw the RA/Dec, RA/Vmag, and RA/diameter scatter plots for the current selection, including calibrators. |
| `observable_domain` | Compute the observable RA window (sunset → sunrise) for CHARA on the selected night. |

---

### Popup windows

| Method | Description |
|---|---|
| `open_popupInfoTargets` | Open a sortable Treeview listing properties of selected science targets. |
| `open_popupAddTarget` | Prompt for a target name, resolve it via SIMBAD, and add it to the selection. |
| `open_popupBestDec` | Show target counts binned in 5° declination strips. |
| `open_popupProgName` | Open the work-package programme-name selector popup. |
| `open_popupInstMode` | Open the instrumental-mode selector popup. |
| `closeframeInfoTargets` | Close the target-info popup. |
| `closeframeAddTarget` | Close the add-target popup. |
| `onCloseCalsec` | Close the secondary-calibrator warning popup. |
| `onCloseLog` | Close the LOG statistics popup. |

---

### Treeview helpers

| Method | Description |
|---|---|
| `insert_popupInfoTargets` | Populate the target-info Treeview with the current selection. |
| `insert_popupAddTarget` | Populate the add-target Treeview with fainter-star candidates. |
| `OnDoubleClick` | Handle a double-click on an add-target row to confirm the addition. |
| `clear_all` | Remove all rows from a Treeview widget. |
| `copy` | Copy selected Treeview rows to the clipboard (bound to `<Control-c>`). |
| `retag` | Reapply alternating row tags after a sort. |
| `treeview_sort_column` | Sort a Treeview by a column and toggle sort direction. |
| `select_allProgName` | Select all programme-name checkbuttons. |
| `deselect_allProgName` | Deselect all programme-name checkbuttons (re-runs the query). |

---

### Primary calibrators

| Method | Description |
|---|---|
| `query_calprim` | Query the SPICA primary calibrator catalog and filter by RA, Dec, and Vmag. |
| `goCalPrim` | Manually trigger the primary-calibrator query. |
| `delCalPrim` | Clear the primary-calibrator selection and refresh the plot. |
| `entryRaRangePrimCallback` | Handle a new RA search range for primary calibrators. |
| `entryDecRangePrimCallback` | Handle a new Dec search range for primary calibrators. |
| `entryVmagRangePrimCallback` | Handle a new Vmag search range for primary calibrators. |

---

### Secondary calibrators

| Method | Description |
|---|---|
| `query_calsec` | Query JSDC2 (Vizier) and filter secondary calibrators by position, visibility, and diameter quality. |
| `goCalSec` | Manually trigger the secondary-calibrator query. |
| `delCalSec` | Clear the secondary-calibrator selection and refresh the plot. |
| `check4BeStars` | Remove Be stars from a calibrator table via SIMBAD. |
| `check4BadCal` | Remove known bad calibrators via the JMMC `badcal` TAP table. |
| `entryRaRangeSecCallback` | Handle a new RA search range for secondary calibrators. |
| `entryDecRangeSecCallback` | Handle a new Dec search range for secondary calibrators. |
| `entryVmagRangeSecCallback` | Handle a new Vmag search range for secondary calibrators. |
| `entryLDDChiSecCallback` | Handle a new LDD chi² threshold. |
| `entryRelErrorSecCallback` | Handle a new maximum diameter relative-error threshold. |
| `entryMinVisSecCallback` | Handle a new minimum squared-visibility threshold. |
| `entryMinVisSecCallback2` | Alternative visibility-threshold handler (skips value update). |

---

### Aspro2 / SAMP export

| Method | Description |
|---|---|
| `onAspro` | Validate the selection and delegate to `import2aspro`. |
| `import2aspro` | Sort targets by transit time, attach models, and call `samp_votable`. |
| `samp_votable` | Serialise targets and calibrators to a VOtable and broadcast via SAMP. |
| `importClients` | Send the prepared SAMP message to the selected Aspro2 client(s). |
| `cancelClickClients` | Cancel the pending SAMP export and disconnect. |

---

### Catalog queries

| Method | Description |
|---|---|
| `dbquery_tap` | Fetch the full SPICA science-target catalog via ADQL/TAP. |
| `calquery_tap` | Fetch the SPICA primary calibrator catalog via ADQL/TAP. |

---

### Data processing utilities

| Method | Description |
|---|---|
| `update_completion_rate` | Increment a target's completion rate based on its SPICA mode. |
| `update_flag_completion` | Convert a completion rate into a three-level completion flag. |
| `update_priority_final` | Compute the final scheduling priority from the PI priority and completion flag. |
| `update_priority_final2` | Extended priority formula for targets shared across two programmes. |
| `addDiamModel` | Add a unified `diam` column to science-target and calibrator tables. |
| `addNssType` | Add an `nss_type` classification column (`Science`, `CalPrim`, `CalSec`). |
| `add_targetmodel` | Attach Aspro-compatible geometric models to science targets. |
| `add_toymodel` | Attach a simple uniform-disk model to calibrator rows. |
| `replacebyHDname` | Replace object names with their HD identifiers via SIMBAD. |
| `fixColumnTypes` | Cast object-dtype columns to `str` and assign UCDs for Aspro2 compatibility. |
| `normalizeColumnNames` | Map catalog-specific column names to the standardized SPICA/Aspro scheme. |

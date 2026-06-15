# SPICA-NSS Tool — Code Documentation

**File:** `nss.py`
**Version:** v2025-09
**Author:** dsalabert
**Created:** December 2023

---

## Table of Contents

1. [Overview](#overview)
2. [Dependencies](#dependencies)
3. [Module-Level Constants](#module-level-constants)
4. [Class: `spica_NSS`](#class-spica_nss)
   - [Class-Level Attributes](#class-level-attributes)
   - [GUI Layout and Initialization (`__init__`)](#gui-layout-and-initialization-__init__)
5. [Methods Reference](#methods-reference)
   - [Database Query Methods](#database-query-methods)
   - [Target Selection Methods](#target-selection-methods)
   - [Calibrator Methods](#calibrator-methods)
   - [Visualization Methods](#visualization-methods)
   - [Aspro2 / SAMP Export Methods](#aspro2--samp-export-methods)
   - [Priority and Completion Methods](#priority-and-completion-methods)
   - [Utility / Helper Methods](#utility--helper-methods)
   - [GUI Event Callbacks](#gui-event-callbacks)
   - [Popup / Dialog Methods](#popup--dialog-methods)
   - [Log and Statistics Methods](#log-and-statistics-methods)
6. [Workflow Summary](#workflow-summary)
7. [Default Filter Parameters](#default-filter-parameters)

---

## Overview

`nss.py` implements the **SPICA-NSS** (Night Science Scheduling) graphical application for observation planning at the CHARA interferometric array. The tool provides an interactive Tkinter GUI to:

- Query the SPICA science target database via a TAP/ADQL service
- Filter targets by work package (program name), instrumental mode, priority, declination, and visual magnitude
- Compute the observable RA domain for any given night at CHARA
- Identify and filter primary and secondary interferometric calibrators
- Visualize sky distribution and magnitude/diameter statistics in real time
- Export observing lists to the JMMC Aspro2 planning tool via the SAMP protocol

The entire application is encapsulated in a single class, `spica_NSS`, which is instantiated at the bottom of the file when the script is run directly.

---

## Dependencies

| Package | Purpose |
|---|---|
| `tkinter` / `ttk` | GUI framework (widgets, dialogs, treeviews) |
| `tkmacosx` | macOS-specific Button widget (imported only on Darwin) |
| `astropy` | Coordinate transforms, VOtable I/O, SAMP client, time handling |
| `astroplan` | Observer site and sunset/sunrise computation |
| `astroquery.simbad` | Querying star names and object types from SIMBAD |
| `astroquery.vizier` | Querying the JSDC2 secondary calibrator catalog |
| `pyvo` | TAP service client for SPICA/JMMC database queries |
| `a2p2.jmmc` | JMMC model/catalog classes for Aspro2-compatible geometric models |
| `numpy` / `numpy.ma` | Numerical operations and masked array handling |
| `matplotlib` | Embedded science plots (RA–Dec, magnitude, diameter) |
| `scipy.special.jv` | Bessel function for squared visibility calculations |
| `json` | Parsing JSON-encoded stellar models from the database |
| `pathlib`, `tempfile` | Temporary VOtable file creation for SAMP transfer |

---

## Module-Level Constants

```python
oidbTapUrl   = "http://tap.jmmc.fr/vollt/tap/"  # OIDB TAP endpoint (currently unused)
tapServerUrl = "http://tap.jmmc.fr/vollt/tap/"  # Active TAP server for SPICA-DB and calibrator queries
nssVersion   = "v2025-09"                        # Application version string
```

---

## Class: `spica_NSS`

```
class spica_NSS
```

Main application class. Contains all GUI construction, database access, filtering logic, visualization, and Aspro2 export functionality.

---

### Class-Level Attributes

These are defined as class-level variables (shared defaults) and can be overridden per instance at runtime through the GUI.

#### State flags

| Attribute | Type | Default | Description |
|---|---|---|---|
| `CalsecOpened` | `bool` | `False` | Whether the secondary calibrator warning popup is open |
| `LogOpened` | `bool` | `False` | Whether the log/statistics window is open |
| `popupInfoTargets` | `bool` | `False` | Whether the target list popup is open |
| `popupAddTarget` | `bool` | `False` | Whether the "add a star" popup is open |
| `popupBestDec` | `bool` | `False` | Whether the best-declination popup is open |
| `topSAMP` | `bool` | `False` | Whether the SAMP/Aspro2 confirmation dialog is open |

#### Catalogs and indices

| Attribute | Type | Default | Description |
|---|---|---|---|
| `date` | `str` | Today's date | Observing date in `YYYY-MM-DD` format |
| `spica_catg` | `Table` or `None` | `None` | Currently loaded SPICA science target catalog |
| `index_AddTarget` | `list` or `None` | `None` | Indices of manually added targets |
| `calprim_catg` | `Table` or `None` | `None` | Primary calibrator catalog |
| `indexList_CalPrim` | `list` or `None` | `None` | Indices of selected primary calibrators |
| `calsec_catg` | `Table` or `None` | `None` | Secondary calibrator catalog |
| `indexList_CalSec` | `ndarray` or `None` | `None` | Indices of selected secondary calibrators |

#### Science target filter defaults

| Attribute | Default | Description |
|---|---|---|
| `decmin` | `-30.0` | Minimum declination (degrees) |
| `decmax` | `90.0` | Maximum declination (degrees) |
| `vmagmin` | `-3.0` | Minimum visual magnitude |
| `vmagmax` | `13.0` | Maximum visual magnitude |

#### Primary calibrator search defaults

| Attribute | Default | Description |
|---|---|---|
| `rarangeprim` | `60` | RA search radius around targets (arcmin) |
| `decrangeprim` | `5` | Dec search radius around median target Dec (degrees) |
| `vmagrangeprim` | `2` | Magnitude search window around median target Vmag |

#### Secondary calibrator search defaults

| Attribute | Default | Description |
|---|---|---|
| `rarangesec` | `60` | RA search radius (arcmin) |
| `decrangesec` | `2` | Dec search radius (degrees) |
| `vmagrangesec` | `2` | Magnitude search window |
| `lddchisec` | `2` | Maximum allowed LDD Chi² for JSDC2 sources |
| `relerrorsec` | `10` | Maximum allowed relative error on limb-darkened diameter (%) |
| `minvissec` | `0.7` | Minimum squared visibility at maximum baseline |
| `maxbaseline` | `330` | Maximum baseline length for visibility computation (meters) |

---

### GUI Layout and Initialization (`__init__`)

```python
def __init__(self)
```

Constructs the main Tkinter window and all its sub-frames, widgets, and layout grid. No parameters are required.

**Initializes:**
- The root `Tk()` window with title, background, and font settings
- All `StringVar` / `IntVar` instances used for widget bindings
- The following labeled frames (rows in the main grid):
  - `FrameDate` — Date entry and display
  - `FrameActions` — Main action buttons (QUERY, BEST_DEC, INFO_TARGETS, SEND2ASPRO, RESET, QUIT, LOG)
  - `FrameWorkPackages` — Program name checkboxes (S01–S08)
  - `FrameInstModes` — Instrument mode checkboxes (DIA, DLD, IMA, TMP, SPI)
  - `FramePriorities` — Final priority checkboxes (1–4)
  - `FrameObjects` — Declination and magnitude range entries with mean display
  - `FrameAddTarget` — "Add a star" button
  - `FrameCalPrims` — Primary calibrator search parameter entries and UNDO button
  - `FrameCalSecs` — Secondary calibrator search parameter entries and baseline slider
  - `FrameLog` — Log/statistics display area
- An embedded `matplotlib` figure with three sub-plots (RA–Dec, Vmag, Diam) attached to the window via `FigureCanvasTkAgg`
- The Tkinter main event loop (`self.root.mainloop()`) at the end of `__init__`

---

## Methods Reference

### Database Query Methods

---

#### `dbquery_tap()`

```python
def dbquery_tap(self) -> astropy.table.Table
```

Queries the full SPICA science target catalog from the JMMC TAP service using ADQL (`SELECT * FROM spica`). Masked values in `completion_rate` and `priority_pi` columns are replaced with numeric defaults to prevent downstream filtering errors.

**Returns:** `astropy.table.Table` — Full SPICA science target catalog.

---

#### `calquery_tap()`

```python
def calquery_tap(self) -> astropy.table.Table
```

Queries the SPICA primary calibrator catalog (`spica_calprim`) from the TAP service.

**Returns:** `astropy.table.Table` — Primary calibrator catalog.

---

### Target Selection Methods

---

#### `onQuery()`

```python
def onQuery(self) -> None
```

Main entry point triggered by the **QUERY_CATALOG** button. Performs the following steps:

1. Closes any open popups (info targets, best-dec)
2. Resets all active selections and index lists
3. Instantiates the CHARA observer via `astroplan`
4. Computes the observable RA domain for the chosen date (sunset → sunrise)
5. On first call: fetches the full SPICA catalog via TAP and caches it; on subsequent calls: reuses the cache and resets filter widgets to defaults
6. Filters targets to those observable at CHARA on the given night using a geometric horizon angle of 70°
7. Activates the relevant program name, instrument mode, and priority checkboxes
8. Replaces star names with HD identifiers where available
9. Triggers `plot_radec()` to refresh the visualization

---

#### `onReset()`

```python
def onReset(self) -> None
```

Resets the full query state by re-initializing `self.iter = 0` and calling `onQuery()`. All filter widgets are returned to their default values.

---

#### `getSelectedTargets()`

```python
def getSelectedTargets(self) -> None
```

Computes the active target list by intersecting all active filter index arrays:

- `indexProgName` — selected program names
- `indexInstMode` — selected instrument modes
- `indexFinalPriority` — selected priority levels
- `indexDecMin`, `indexDecMax` — declination bounds
- `indexVmagMin`, `indexVmagMax` — magnitude bounds

Stores the result in `self.indexList_Targets`. Manually added targets (`self.index_AddTarget`) are appended after the intersection.

---

#### `getAddTarget()`

```python
def getAddTarget(self) -> None
```

Similar to `getSelectedTargets()`, but uses a separate magnitude filter (`indexVmagToAddTarget`) to identify fainter candidate stars that could be manually added to the selection.

---

#### `plotSelectedProgName()`

```python
def plotSelectedProgName(self) -> None
```

Callback for program name checkboxes. Updates `self.indexProgName` based on checked programs (searching both `progname` and `progname2` columns), then refreshes the target list and plot.

---

#### `plotSelectedInstMode()`

```python
def plotSelectedInstMode(self) -> None
```

Callback for instrument mode checkboxes. Filters `self.indexInstMode` by matching the `spica_mode` column, then refreshes the target list and plot.

---

#### `plotSelectedFinalPriority()`

```python
def plotSelectedFinalPriority(self) -> None
```

Callback for priority checkboxes. Filters `self.indexFinalPriority` by matching the `priority_final` column, then refreshes the target list and plot.

---

### Calibrator Methods

---

#### `query_calprim()`

```python
def query_calprim(self) -> None
```

Searches for primary interferometric calibrators for the currently selected science targets.

**Procedure:**
1. Computes RA, Dec, and Vmag search windows centred on the selected target sample
2. Queries `spica_calprim` via TAP (`calquery_tap()`)
3. Filters by nighttime observability (RA between sunset and sunrise)
4. Filters by Dec and Vmag ranges
5. Removes any calibrators that coincide with science targets (within 5 arcsec)
6. Replaces names with HD identifiers
7. Stores results in `self.calprim_catg` and `self.indexList_CalPrim`
8. Refreshes the plot

---

#### `query_calsec()`

```python
def query_calsec(self) -> None
```

Searches for secondary interferometric calibrators using the JSDC2 catalog (VizieR `II/346/jsdc_v2`).

**Procedure:**
1. Computes RA, Dec, and Vmag search windows
2. Queries JSDC2 with filters on Dec, Vmag, `CalFlag = 0`, and LDD Chi² < `lddchisec`
3. Filters by nighttime observability
4. Computes squared visibility `vis²` at the maximum baseline using the Bessel function formula:
   ```
   vis² = (2·J₁(z) / z)²   where z = π·d·B / λ
   ```
5. Filters by relative diameter error and minimum visibility
6. Removes stars already in the science target list or in the primary calibrator list (5 arcsec crossmatch)
7. Removes Be stars (`check4BeStars`) and known bad calibrators (`check4BadCal`)
8. Warns the user if more than 50 secondary calibrators remain (cap enforced at 50)
9. Stores results in `self.calsec_catg` and `self.indexList_CalSec`

---

#### `check4BeStars(tableObjects)`

```python
def check4BeStars(self, tableObjects: Table) -> Table
```

Queries SIMBAD for each star in the table. Removes any object classified as `Be*` (Be star), which are unsuitable as interferometric calibrators.

**Parameters:**
- `tableObjects` — Astropy Table with a `Name` column

**Returns:** Table with Be stars removed.

---

#### `check4BadCal(tableObjects)`

```python
def check4BadCal(self, tableObjects: Table) -> Table
```

Queries the `badcal` table from the TAP service and cross-matches it with the input catalog using a 5 arcsec sky matching radius. Removes any objects flagged as bad calibrators.

**Parameters:**
- `tableObjects` — Astropy Table with `ra` and `dec` columns

**Returns:** Table with known bad calibrators removed.

---

#### `delCalPrim()`

```python
def delCalPrim(self) -> None
```

Clears the primary calibrator selection and refreshes the plot. Bound to the **UNDO** button in the primary calibrators frame.

---

#### `delCalSec()`

```python
def delCalSec(self) -> None
```

Clears the secondary calibrator selection and refreshes the plot.

---

### Visualization Methods

---

#### `plot_radec()`

```python
def plot_radec(self) -> None
```

Redraws the embedded matplotlib figure, which consists of three vertically stacked subplots:

- **Top (ax[0]):** LST vs. Declination scatter plot. Point size scales inversely with Vmag. Colors encode final priority (1=red, 2=blue, 3=green, 4=orange). Shape encodes instrument mode. Vertical dashed lines mark sunset and sunrise LST.
- **Middle (ax[1]):** LST vs. visual magnitude (Vmag).
- **Bottom (ax[2]):** LST vs. angular diameter (Diam, in mas).

Science targets are shown by priority level; primary calibrators appear in blue; secondary calibrators appear in cyan. A legend box summarizes the count of targets per priority per mode (DIA, DLD, IMA, TMP, SPI).

The x-axis is the Local Sidereal Time (LST) in hours, adjusted so that midnight appears at the centre.

---

#### `observable_domain()`

```python
def observable_domain(self) -> tuple[float, float]
```

Computes the observable right ascension range for the CHARA observatory on the configured date.

**Procedure:**
1. Retrieves sunset and sunrise UTC times using `astroplan`
2. Converts both to Local Mean Sidereal Time (LMST)
3. Converts LMST to degrees (×15)

**Returns:** `(alpha_sun_set, alpha_sun_rise)` — RA boundaries in degrees.

---

### Aspro2 / SAMP Export Methods

---

#### `onAspro()`

```python
def onAspro(self) -> None
```

Triggered by the **SEND2ASPRO** button. Calls `samp_votable()` with the selected science targets, primary calibrators, and secondary calibrators to initiate the Aspro2 export workflow.

---

#### `samp_votable(targets, calibrators1=None, calibrators2=None)`

```python
def samp_votable(
    self,
    targets: Table,
    calibrators1: Table = None,
    calibrators2: Table = None
) -> None
```

Converts the selected targets and calibrators into a VOtable and transmits it to a running Aspro2 instance via SAMP.

**Procedure:**
1. Connects as a SAMP client named `"SPICA-NSS"`
2. Normalizes column names for all tables (`normalizeColumnNames`)
3. Assigns group labels (`calprim`, `calsecond`, `priority_pi=N`) for Aspro2 color coding
4. Enriches comments with PI name, program name, and mode
5. Attaches geometric stellar models (`add_targetmodel`, `add_toymodel`)
6. Stacks science targets, primary calibrators, and secondary calibrators into a single table
7. Sorts the merged list by meridian transit time at CHARA
8. Writes the VOtable to a temporary file with observation metadata (interferometer, period, instrument, date, baseline configuration)
9. Detects connected Aspro2 SAMP clients and shows a confirmation dialog:
   - 0 clients → error message
   - 1 client → simple Yes/No confirmation
   - Multiple clients → checkbox selector per client ID

---

#### `importClients()`

```python
def importClients(self) -> None
```

Sends the prepared SAMP message (VOtable URI) to the selected Aspro2 client(s) and disconnects.

---

#### `cancelClickClients()`

```python
def cancelClickClients(self) -> None
```

Cancels the SAMP export and disconnects the SAMP client without sending data.

---

#### `add_targetmodel(data)`

```python
def add_targetmodel(self, data: Table) -> Table
```

Attaches Aspro2-compatible XML geometric models to science targets. If the target has a JSON-encoded model in the database, it is converted to Aspro's internal `_model` format. If no model is present, a simple uniform-disk model is created from the angular diameter column `diam`.

---

#### `add_toymodel(data)`

```python
def add_toymodel(self, data: Table) -> Table
```

Creates a simple disk model for calibrators using the `ld_jsdc2` (limb-darkened diameter from JSDC2) column and attaches it to the table in the Aspro2 `model` column format.

---

### Priority and Completion Methods

---

#### `update_completion_rate(completion_rate, spica_mode)`

```python
def update_completion_rate(self, completion_rate: float, spica_mode: str) -> float
```

Increments the completion rate by a mode-dependent step to simulate one additional observation (used in the DRS simulation):

| Mode | Increment |
|---|---|
| DIA | +0.50 |
| DLD | +0.25 |
| IMA / TMP | +0.10 |
| SPI | +0.20 |

**Returns:** Updated completion rate (capped at the existing completion logic).

---

#### `update_flag_completion(completion_rate, spica_mode)`

```python
def update_flag_completion(self, completion_rate: float, spica_mode: str) -> int
```

Converts a completion rate into a discrete completion flag used to determine scheduling priority:

| Flag | Meaning |
|---|---|
| `1` | Sufficiently completed |
| `2` | Partially completed |
| `3` | Not started or SPI mode |

The thresholds depend on the mode (e.g., DIA: flag=1 if completion_rate > 0; DLD/IMA/TMP: flag=1 if ≥ 0.5, flag=2 if 0 < rate < 0.5).

---

#### `update_priority_final(flag_completion, priority_pi)`

```python
def update_priority_final(self, flag_completion: int, priority_pi: int) -> int
```

Computes a final priority integer (1–4) by combining the PI-assigned priority with the completion flag. Lower numbers indicate higher scheduling priority.

---

#### `update_priority_final2(flag_completion, priority_pi, progname2)`

```python
def update_priority_final2(self, flag_completion: int, priority_pi: int, progname2) -> int
```

Extended version of `update_priority_final` that also accounts for shared targets appearing in a second program (`progname2`). Used when a target belongs to two SPICA work packages simultaneously.

---

### Utility / Helper Methods

---

#### `replacebyHDname(tableObjects)`

```python
def replacebyHDname(self, tableObjects: Table) -> Table
```

For each object in the input catalog, queries SIMBAD for alternative identifiers. If an HD designation is found and the current name is not already an HD identifier, the name is replaced. Works on columns named `target_main_id`, `name`, or `Name`.

**Returns:** Table with HD-prefixed names where available.

---

#### `fixColumnTypes(ts)`

```python
def fixColumnTypes(self, ts: list[Table]) -> None
```

Converts `object`-dtype columns to `str` and assigns IVOA UCD metadata strings to standard columns (`ra`, `dec`, `vmag`, `hmag`, etc.). This ensures compatibility when stacking tables from different catalog sources and when writing VOtables for Aspro2.

---

#### `normalizeColumnNames(input, colNames=...)`

```python
def normalizeColumnNames(self, input: Table, colNames: dict = {...}) -> Table
```

Renames catalog-specific column names to a standardized SPICA/Aspro naming scheme. Allows tables from SPICA-DB, JSDC2, and calibrator catalogs to be merged consistently. Applies `fixColumnTypes` to the output.

**Default mapping includes:** `ra`, `dec`, `name→target_main_id`, `ldd→ld_jsdc2`, `vmag`, `hmag`, `diam`, `sptype→spt`, etc.

---

#### `addDiam(targets, calibrators1=None, calibrators2=None)`

```python
def addDiam(self, targets, calibrators1=None, calibrators2=None)
```

Adds a unified `diam` column to science targets (extracted from the JSON model field, supporting `disk` and `elong_disk` types), primary calibrators (from `ldd` column), and secondary calibrators (from `UDDR` column). Returns all three tables.

---

#### `addNssType(targets, calibrators1=None, calibrators2=None)`

```python
def addNssType(self, targets, calibrators1=None, calibrators2=None)
```

Appends an `nss_type` column to each table to distinguish object roles in the exported catalog: `"Science"`, `"CalPrim"`, or `"CalSec"`.

---

#### `retag(theTreeToSort)`

```python
def retag(self, theTreeToSort: ttk.Treeview) -> None
```

Reapplies alternating `oddrow`/`evenrow` background tags to a Treeview after sorting, preserving the striped row appearance.

---

#### `treeview_sort_column(tv, col, reverse)`

```python
def treeview_sort_column(self, tv: ttk.Treeview, col: str, reverse: bool) -> None
```

Sorts a Treeview column in ascending or descending order when the column header is clicked. Integer sorting is applied to the `Count` column; string sorting is used for all others.

---

#### `clear_all(tree)`

```python
def clear_all(self, tree: ttk.Treeview) -> None
```

Removes all rows from the given Treeview widget.

---

#### `copy(event)`

```python
def copy(self, event) -> None
```

Copies selected Treeview row(s) to the system clipboard as tab-separated values, triggered by `Ctrl+C`.

---

#### `onDRS()`

```python
def onDRS(self) -> None
```

Simulates a Data Reduction System (DRS) pass: increments the `completion_rate` and `ob_refs` of selected targets and resets `qcs_flag` to `0`. Updates the in-memory catalog cache (`spica_catg_jmmc`) accordingly. This is a development/testing feature.

---

### GUI Event Callbacks

These short methods are bound to GUI entry widgets and fire `<Return>` key events or slider releases. They update the corresponding instance attribute and trigger re-queries or re-plots.

| Method | Widget | Action |
|---|---|---|
| `entryDateCallback(strDate)` | Date entry | Sets `self.date`, updates label, calls `onQuery()` |
| `entryDecMinCallback(strDecMin)` | DEC_MIN entry | Updates `self.decmin`, refreshes index and plot |
| `entryDecMaxCallback(strDecMax)` | DEC_MAX entry | Updates `self.decmax`, refreshes index and plot |
| `entryVmagMinCallback(strVmagMin)` | VMAG_MIN entry | Updates `self.vmagmin`, refreshes index and plot |
| `entryVmagMaxCallback(strVmagMax)` | VMAG_MAX entry | Updates `self.vmagmax`, refreshes index and plot |
| `entryRaRangePrimCallback(...)` | RA range (prim) entry | Updates `self.rarangeprim`, calls `query_calprim()` |
| `entryDecRangePrimCallback(...)` | DEC range (prim) entry | Updates `self.decrangeprim`, calls `query_calprim()` |
| `entryVmagRangePrimCallback(...)` | Vmag range (prim) entry | Updates `self.vmagrangeprim`, calls `query_calprim()` |
| `entryRaRangeSecCallback(...)` | RA range (sec) entry | Updates `self.rarangesec`, calls `query_calsec()` |
| `entryDecRangeSecCallback(...)` | DEC range (sec) entry | Updates `self.decrangesec`, calls `query_calsec()` |
| `entryVmagRangeSecCallback(...)` | Vmag range (sec) entry | Updates `self.vmagrangesec`, calls `query_calsec()` |
| `entryLDDChiSecCallback(...)` | Max LDD Chi² entry | Updates `self.lddchisec`, calls `query_calsec()` |
| `entryRelErrorSecCallback(...)` | Max rel. error entry | Updates `self.relerrorsec`, calls `query_calsec()` |
| `entryMinVisSecCallback(...)` | Min vis² entry | Updates `self.minvissec`, calls `query_calsec()` |
| `entryMinVisSecCallback2(...)` | Max baseline slider | Triggers `query_calsec()` on slider release |
| `goCalPrim()` | (internal) | Alias: calls `query_calprim()` |
| `goCalSec()` | (internal) | Alias: calls `query_calsec()` |

---

### Popup / Dialog Methods

---

#### `open_popupAddTarget()`

```python
def open_popupAddTarget(self) -> None
```

Prompts the user to enter a star name, queries SIMBAD to resolve it, and cross-matches against the loaded SPICA catalog. If the target is found in the catalog but not yet selected, the user is asked to confirm its addition. If multiple entries match (same star in different programs), a Treeview dialog lets the user double-click to pick the desired entry.

---

#### `OnDoubleClick(event)`

```python
def OnDoubleClick(self, event) -> None
```

Handles double-click events on the add-target Treeview. Confirms the selection, appends the index to `self.index_AddTarget`, and refreshes the plot and info popup.

---

#### `open_popupInfoTargets()`

```python
def open_popupInfoTargets(self) -> None
```

Opens a sortable Treeview popup listing all currently selected science targets with columns for: SPICA-DB ID, target name, spectral type, program name, mode, priority, completion rate, RA, Dec, diameter (mas), Vmag, and Hmag.

---

#### `insert_popupInfoTargets()`

```python
def insert_popupInfoTargets(self) -> None
```

Populates the info targets Treeview with data from `self.spica_catg[self.indexList_Targets]`, sorted by Vmag. Handles dual-program targets (displaying combined progname and priority strings).

---

#### `insert_popupAddTarget()`

```python
def insert_popupAddTarget(self) -> None
```

Populates the add-target Treeview with candidate faint stars from `self.indexList_AddTarget`, sorted by Vmag.

---

#### `open_popupBestDec()`

```python
def open_popupBestDec(self) -> None
```

Opens a Treeview showing the number of targets in each 5° declination bin between −30° and +90°, within the current Vmag filter. The Count column is sortable. Helps observers choose the best declination range to observe.

---

#### `open_popupProgName()`

```python
def open_popupProgName(self) -> None
```

Opens a secondary popup listing all available program names as checkboxes, with SELECT ALL and RESET controls. Allows finer filtering than the top-level inline checkboxes.

---

#### `open_popupInstMode()`

```python
def open_popupInstMode(self) -> None
```

Opens a secondary popup for instrumental mode selection, dynamically populated from the unique modes present in the current program name selection. SPI mode is always listed last.

---

#### `select_allProgName()` / `deselect_allProgName()`

```python
def select_allProgName(self) -> None
def deselect_allProgName(self) -> None
```

Selects all program name checkboxes, or resets them by re-running `onQuery()`.

---

#### `onCloseCalsec()`

```python
def onCloseCalsec(self) -> None
```

Closes the secondary calibrator overflow warning popup and resets `self.CalsecOpened`.

---

### Log and Statistics Methods

---

#### `onLog()`

```python
def onLog(self) -> None
```

Opens or refreshes the log/statistics window. Displays survey completion statistics broken down by priority level and program, based on the full (unfiltered) SPICA catalog. Queries the JMMC database.

---

#### `onQuit()`

```python
def onQuit(self) -> None
```

Prompts the user for confirmation and quits the application via `self.root.destroy()`.

---

## Workflow Summary

A typical observing preparation session with SPICA-NSS follows this sequence:

```
1. Launch → __init__ builds the GUI
2. Enter a date (or use today's default) → entryDateCallback → onQuery
   └─ Fetches SPICA-DB via TAP
   └─ Computes CHARA nighttime window
   └─ Filters targets observable that night
   └─ Plots RA–Dec / Vmag / Diam charts
3. Use checkboxes to filter by program, mode, priority
4. Adjust DEC and VMAG range sliders as needed
5. Optionally add a custom star with "Add a star" → open_popupAddTarget
6. Click BEST_DEC → open_popupBestDec (declination distribution table)
7. Click INFO_TARGETS → open_popupInfoTargets (sortable target list)
8. Enter RA/Dec/Vmag ranges → query_calprim (primary calibrators appear on plot)
9. Enter search ranges + vis²/error filters → query_calsec (secondary calibrators appear on plot)
10. Click SEND2ASPRO → samp_votable → VOtable broadcast to Aspro2 via SAMP
11. Click QUIT to exit
```

---

## Default Filter Parameters

| Parameter | Default | Unit | Adjustable via |
|---|---|---|---|
| Declination min | −30 | degrees | DEC_MIN entry |
| Declination max | +90 | degrees | DEC_MAX entry |
| Vmag min | −3 | mag | VMAG_MIN entry |
| Vmag max | +13 | mag | VMAG_MAX entry |
| Primary RA range | 60 | arcmin | RA range (prim) entry |
| Primary Dec range | 5 | degrees | DEC range (prim) entry |
| Primary Vmag range | ±2 | mag | Vmag range (prim) entry |
| Secondary RA range | 60 | arcmin | RA range (sec) entry |
| Secondary Dec range | 2 | degrees | DEC range (sec) entry |
| Secondary Vmag range | ±2 | mag | Vmag range (sec) entry |
| Max LDD Chi² | 2 | — | LDD Chi² entry |
| Max relative diameter error | 10 | % | Rel. error entry |
| Min squared visibility | 0.7 | — | Min vis² entry |
| Max baseline | 330 | meters | Baseline slider |

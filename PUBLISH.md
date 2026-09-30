# PUBLISH.md — how to put this deposit on GitHub and archive it on Zenodo

Two parts:

* **Part 1** — create the GitHub repository and push, in Windows PowerShell.
* **Part 2** — link Zenodo, enable the repository, reserve the DOI, and cut the release that mints it.

Target repository: `https://github.com/PlatypusConservation/KI-platypus`
Target release tag: `v1.0.0`

The organisation `PlatypusConservation` exists. The repository does **not** exist yet — Part 1
creates it. Nothing below has been run for you.

Before you start, read **§2.3** if the DOI has to go into the manuscript *before* the release is
live. That changes the order of the steps, and it is easier to do in the right order than to fix
afterwards.

---

# Part 1 — GitHub, in Windows PowerShell

## 1.0 Open PowerShell in the deposit folder

**Do this first, before any git command.** A fresh PowerShell window often opens in
`C:\Windows\System32`. Running `git init` there creates a repository inside a Windows system
folder. It fails harmlessly at the commit step, but it leaves a stray `.git` behind that you then
have to remove.

If that has already happened, remove it and do **not** run the `safe.directory` line git suggests —
that would permanently whitelist a system folder:

```powershell
Remove-Item -LiteralPath "C:\Windows\System32\.git" -Recurse -Force
Test-Path "C:\Windows\System32\.git"      # must print False
```

**Work outside OneDrive.** OneDrive syncs and locks files inside `.git`, which causes intermittent
git errors. Copy the deposit to a plain local folder and work there:

```powershell
$src = "C:\Users\webre\OneDrive - UNSW\E\Platypus\Kangaroo Island\Genetics\Revision\deposit"
$dst = "C:\Users\webre\KI-platypus"
Copy-Item -LiteralPath $src -Destination $dst -Recurse
Set-Location -LiteralPath $dst
```

Confirm you are in the right place before going further. This must list `README.md`, `LICENSE`,
`CITATION.cff`, `PUBLISH.md`, and the `code`, `data`, `figures` and `results` folders:

```powershell
Get-Location
Get-ChildItem -Force | Select-Object Name
```

PowerShell hides dot-files by default; `-Force` shows `.gitignore` and `.zenodo.json`.

If the listing is wrong, stop — you are in the wrong folder, and every command below would act on it.

An elevated ("Run as administrator") window is not needed for any of this, and is what lands you in
System32 in the first place. A normal PowerShell window is the better choice.

## 1.1 Check what is installed

Run these four lines. Each prints something or prints nothing; nothing means "not installed".

```powershell
Get-Command git -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Source
git --version
Get-Command gh  -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Source
gh --version
```

* **`git --version` printed a version** → good, continue.
* **`git` printed nothing** → install Git for Windows from <https://git-scm.com/download/win>, close
  PowerShell, open a new one, `Set-Location` back to the deposit folder, and re-run this section.
  Nothing further works without git.
* **`gh --version` printed a version** → use **Variant A** in §1.5. This is the short path.
* **`gh` printed nothing** → use **Variant B** in §1.6. You will create the repository in a browser.
  Do not install `gh` on account of this; Variant B is only three commands longer.

Now check that git knows who you are. These two must both print a value:

```powershell
git config --global user.name
git config --global user.email
```

If either is empty, set them — the commit will be rejected or attributed to nobody otherwise:

```powershell
git config --global user.name  "Gilad Bino"
git config --global user.email "gilad.bino@unsw.edu.au"
```

Use the email address that is on your GitHub account, or GitHub will not link the commit to you.

## 1.2 Check the folder is a clean starting point

This must **not** already be a git repository. Check:

```powershell
Test-Path ".git"
```

* `False` → good, continue to §1.3.
* `True` → the folder is already a repository. Stop and run `git remote -v` and `git log --oneline -5`
  to see what it is before doing anything else. Do not run `git init` again.

## 1.3 Confirm what will be committed

The `.gitignore` in this folder excludes Python caches, virtual environments and editor droppings.
Everything else goes up, including the 1.8 MB genotype matrix and the figure PDFs. Check the total
size first — GitHub warns above 1 GB per repository and rejects any single file above 100 MB, and
this deposit is far below both:

```powershell
"{0:N1} MB" -f ((Get-ChildItem -Recurse -File -Force |
  Where-Object { $_.FullName -notmatch '\\\.git\\' } |
  Measure-Object -Property Length -Sum).Sum / 1MB)

Get-ChildItem -Recurse -File -Force |
  Where-Object { $_.Length -gt 20MB } |
  Select-Object FullName, @{n='MB';e={[math]::Round($_.Length/1MB,1)}}
```

The first command should print roughly **11 MB**. The second should print nothing. If it lists a
file, stop and decide whether that file belongs in the deposit before you commit it.

## 1.4 Initialise, stage, commit

```powershell
git init -b main
git add -A
git status
```

Read the `git status` output before committing. It should list `README.md`, `LICENSE`,
`CITATION.cff`, `PUBLISH.md`, `.zenodo.json`, `.gitignore`, and the files under `code/`, `data/`,
`figures/` and `results/` — and nothing else. Then:

```powershell
git commit -m "Code and data for the Kangaroo Island platypus genomic diversity paper (v1.0.0)"
git log --oneline
```

If `git init -b main` failed with an error about `-b` (git older than 2.28), do this instead:

```powershell
git init
git add -A
git commit -m "Code and data for the Kangaroo Island platypus genomic diversity paper (v1.0.0)"
git branch -M main
```

Confirm the branch is named `main` — GitHub's default, and what the rest of this file assumes:

```powershell
git branch --show-current
```

## 1.5 Variant A — with the GitHub CLI (`gh`)

Use this only if `gh --version` printed a version in §1.1.

First authenticate. Check whether you already are:

```powershell
gh auth status
```

If it says you are not logged in, run the interactive login and follow its prompts (choose
`GitHub.com`, then `HTTPS`, then authenticate in the browser):

```powershell
gh auth login
```

Confirm you can write to the organisation. This must list `PlatypusConservation`:

```powershell
gh org list
```

If it does not, you are signed in as an account without access to that organisation, or the
authentication is missing the `read:org` scope — run `gh auth refresh -s read:org,repo` and re-check.

Now create the repository and push in one step:

```powershell
gh repo create PlatypusConservation/KI-platypus --public --source . --remote origin --push --description "Code and data for Bino et al., genome-wide diversity in the isolated, introduced Kangaroo Island platypus (Ornithorhynchus anatinus)"
```

Verify:

```powershell
git remote -v
git status
gh repo view PlatypusConservation/KI-platypus --web
```

`git status` should say `Your branch is up to date with 'origin/main'`. The browser should show your
files and the rendered `README.md`.

**Then stop here and go to Part 2 before cutting the release.** Zenodo has to be watching the
repository *before* the release exists, or the release is not archived.

When Part 2 §2.2 tells you to come back, cut the release:

```powershell
gh release create v1.0.0 --title "v1.0.0 — Conservation Genetics submission" --notes "Code and data accompanying Bino, Hawke, Baring and Gongora, 'Drifting alone: genome-wide diversity in the isolated, introduced Kangaroo Island platypus (Ornithorhynchus anatinus)', Conservation Genetics. Analysis dataset: 4,002 biallelic SNPs, sex-linked loci excluded, scored in all 222 individuals with no missing genotypes. See README.md for the manifest, the reproduction instructions and the data-provenance note."
gh release view v1.0.0 --web
```

## 1.6 Variant B — plain `git`, repository created in the browser

Use this if `gh` is not installed.

**Step 1 — create the empty repository on GitHub.**

1. Go to <https://github.com/organizations/PlatypusConservation/repositories/new>.
2. **Owner**: `PlatypusConservation`. **Repository name**: `KI-platypus`.
3. **Description**: `Code and data for Bino et al., genome-wide diversity in the isolated, introduced Kangaroo Island platypus (Ornithorhynchus anatinus)`
4. Select **Public**. Zenodo cannot archive a private repository.
5. Under "Initialize this repository with", leave **every box unticked** — no README, no .gitignore,
   no licence. This matters: if GitHub creates a first commit, your push in step 2 is rejected as a
   non-fast-forward and you have to merge unrelated histories to recover.
6. Click **Create repository**.

GitHub then shows a page headed "…or push an existing repository from the command line". Ignore its
commands and use the ones below, which match what you have already done.

**Step 2 — add the remote and push.**

```powershell
git remote add origin https://github.com/PlatypusConservation/KI-platypus.git
git remote -v
git push -u origin main
```

If you are asked to sign in, a browser window opens (Git Credential Manager) — sign in there. If you
are asked for a password in the terminal instead, GitHub does not accept your account password: use a
personal access token as the password, created at
<https://github.com/settings/tokens> with the `repo` scope.

Verify:

```powershell
git status
Start-Process "https://github.com/PlatypusConservation/KI-platypus"
```

`git status` should say `Your branch is up to date with 'origin/main'`.

**Then stop here and go to Part 2 before cutting the release.** Zenodo has to be watching the
repository *before* the release exists.

**Step 3 — when Part 2 §2.2 sends you back, cut the tagged release.**

Create and push an annotated tag:

```powershell
git tag -a v1.0.0 -m "v1.0.0 — Conservation Genetics submission"
git push origin v1.0.0
git tag --list
```

**A pushed tag alone does not trigger Zenodo.** Zenodo listens for a *published GitHub Release*.
Turn the tag into one:

1. Go to <https://github.com/PlatypusConservation/KI-platypus/releases/new>.
2. **Choose a tag**: pick the existing `v1.0.0` from the dropdown. Do not type a new one.
3. **Release title**: `v1.0.0 — Conservation Genetics submission`
4. **Describe this release**: paste

   > Code and data accompanying Bino, Hawke, Baring and Gongora, "Drifting alone: genome-wide
   > diversity in the isolated, introduced Kangaroo Island platypus (*Ornithorhynchus anatinus*)",
   > *Conservation Genetics*. Analysis dataset: 4,002 biallelic SNPs, sex-linked loci excluded,
   > scored in all 222 individuals with no missing genotypes. See README.md for the manifest, the
   > reproduction instructions and the data-provenance note.

5. Leave **Set as a pre-release** unticked. Zenodo ignores pre-releases by default.
6. Click **Publish release**.

## 1.7 Later changes

Any further edit — for instance writing the real DOI into `README.md` and `CITATION.cff` — goes up
the same way, and needs a new release if you want Zenodo to archive it:

```powershell
git add -A
git commit -m "Insert the Zenodo DOI in README.md and CITATION.cff"
git push
```

Then cut `v1.0.1` exactly as you cut `v1.0.0` (§1.5 or §1.6 step 3, with the tag changed).

---

# Part 2 — Zenodo: archive the repository and mint the DOI

## 2.1 Link your GitHub account to Zenodo

1. Go to <https://zenodo.org> and click **Sign up** (or **Log in**).
2. Choose **Sign up / Log in with GitHub**. Using GitHub as the sign-in is the simplest route; if you
   already have a Zenodo account with a different login, log into that one instead and connect GitHub
   at <https://zenodo.org/account/settings/linkedaccounts/>.
3. Authorise Zenodo when GitHub asks. Zenodo requests permission to read your repositories and to
   install a release webhook. To see repositories owned by `PlatypusConservation` you must also
   **grant Zenodo access to that organisation** on the GitHub authorisation screen — if the
   organisation shows a **Grant** button next to it, click it. Without that grant the repository will
   not appear in the next step.

## 2.2 Switch the repository on — before the release exists

1. Go to <https://zenodo.org/account/settings/github/>.
2. You will see a list of repositories you can archive, grouped by owner. Find
   **`PlatypusConservation/KI-platypus`** and set its toggle to **ON**.
3. If it is not listed, click **Sync now** at the top right of that page. If it is still missing, the
   organisation grant in §2.1 step 3 did not happen — go to
   <https://github.com/settings/applications>, open **Zenodo**, and grant access to
   `PlatypusConservation`, then **Sync now** again.

Flipping the toggle installs a webhook on the repository. **Only releases published after the toggle
is ON are archived.** A release you published earlier is invisible to Zenodo and cannot be picked up
retrospectively; you would have to cut another one.

**Now go back to Part 1 §1.5 or §1.6 step 3 and publish the `v1.0.0` release.**

Within a minute or two of publishing it, a new record appears in your Zenodo uploads. The record's
files are a zip snapshot of the repository at the tag, and its metadata is read from `.zenodo.json`
in this folder — which is why that file exists and why it is worth getting right before the release.
Open the record and check the title, the four creators, the description and the licence. Add each
creator's **ORCID** here: `.zenodo.json` deliberately omits the `orcid` fields, because an empty ORCID
fails Zenodo's validation, so they have to be typed in once. Click **Publish** / **Save** when done.

## 2.3 Reserving the DOI *before* the release goes live

You asked for the real DOI to be in the manuscript before anything is public. Zenodo supports that,
but not on a record the GitHub integration creates for you — that record only exists once the release
is published. So there are two routes. Pick one **before** you cut the release.

### Route A — reserve the DOI on a manual Zenodo upload (use this if the DOI must go in the manuscript now)

1. On Zenodo, click **New upload** (<https://zenodo.org/uploads/new>).
2. In the upload form, the first field block is **Digital Object Identifier**. Leave the answer to
   "Do you already have a DOI for this upload?" as **No** — that is the default — and click the
   **Reserve DOI** button that sits under that field. (Older Zenodo labels the same button
   **Get a DOI now!**.)
3. Zenodo immediately displays the reserved DOI, in the form `10.5281/zenodo.NNNNNNN`, in that same
   field. It is yours from that moment and will not be given to anyone else. The record stays an
   unpublished draft and nothing is public.
4. Fill in the rest of the form from `.zenodo.json`: **Resource type** = *Dataset*; the title;
   the four creators with affiliations and ORCIDs; the description; **Licence** = *Creative Commons
   Zero v1.0 Universal*; the keywords; and under **Related works** add
   `https://github.com/PlatypusConservation/KI-platypus` as *is supplement to* and
   `10.5281/zenodo.7039778` as *references*.
5. Upload the files. Make a zip that contains the deposit and not the git metadata, by naming the
   items explicitly — run this from the deposit folder:

   ```powershell
   Compress-Archive -Path code, data, figures, results, README.md, LICENSE, CITATION.cff, PUBLISH.md, .zenodo.json, .gitignore -DestinationPath "$HOME\Downloads\KI-platypus-v1.0.0.zip" -Force
   "{0:N1} MB" -f ((Get-Item "$HOME\Downloads\KI-platypus-v1.0.0.zip").Length / 1MB)
   ```

   Then drag that zip into the Zenodo upload form. (If you have already pushed and tagged, you can
   instead download GitHub's own snapshot of the tag, which excludes `.git` for you:
   `https://github.com/PlatypusConservation/KI-platypus/archive/refs/tags/v1.0.0.zip`.)
6. Click **Save draft**, not **Publish**. Put the reserved DOI into the manuscript, into `README.md`
   and into `CITATION.cff` (three placeholders, all reading `10.5281/zenodo.XXXXXXX`).
7. When the paper is accepted, come back and click **Publish**. The reserved DOI activates and
   resolves. Until then the DOI is registered but does not resolve, which is normal and expected for
   a reserved DOI.
8. The GitHub integration can still be used afterwards, for the code's own home and for later
   versions; it will create a *separate* Zenodo record with its own DOI. If you would rather have one
   record only, use Route A alone and do not switch the repository on in §2.2 — you keep the GitHub
   repository as the working home of the code and Zenodo as the citable archive.

### Route B — let the GitHub release mint the DOI (simpler, if the DOI can be added at proof stage)

Cut `v1.0.0` as in Part 1, let Zenodo create the record, take the DOI it mints, and put it into the
manuscript at the proof stage. Then write the DOI into `README.md` and `CITATION.cff` and cut
`v1.0.1` (§1.7) so the repository itself carries its own DOI. The concept DOI (§2.4) is stable across
both releases, so the DOI you gave the journal stays correct.

Route B is less work and is what most authors do. Route A is the one that answers "the DOI must be in
the manuscript before anything is public".

## 2.4 Concept DOI versus version DOI — and which to cite

Zenodo mints **two** DOIs, and they differ by one digit, so it is easy to use the wrong one.

| | What it identifies | Behaviour |
|---|---|---|
| **Concept DOI** (also "all versions") | The record as a work, across every version | Always resolves to the **newest** version. Never changes. |
| **Version DOI** | One specific release, e.g. `v1.0.0` | Always resolves to **that** release, frozen. A new release gets a new one. |

On the Zenodo record page, look at the right-hand **Versions** panel. The version DOI is at the top,
beside the version you are looking at; the concept DOI is the one labelled **"Cite all versions?"**
or **"All versions"** just below the version list.

**Cite the concept DOI in the paper.** Reasons: if you fix a typo in the README, or add the DOI to
`CITATION.cff` and cut `v1.0.1`, a reader following the concept DOI lands on the corrected deposit,
while a reader following the `v1.0.0` version DOI lands on the version that lacks the fix. The
concept DOI is also what stays correct if you later add the peer-review-driven version of the code.
Name the version in the reference text so the reader knows what you used, and let the DOI carry them
to the latest:

> Bino, G., Hawke, T., Baring, R. & Gongora, J. (2026). *Code and data for "Drifting alone:
> genome-wide diversity in the isolated, introduced Kangaroo Island platypus (Ornithorhynchus
> anatinus)"* (v1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.NNNNNNN

Use the **version** DOI only where the exact bytes matter — a reproducibility appendix, or a
reviewer's request to see precisely what was run.

Under Route A there is only one version, so the concept and version DOIs point at the same files;
citing the concept DOI still future-proofs the reference.

## 2.5 Finish up

1. Put the DOI into the manuscript's **Data availability** and **Code availability** statements. A
   sentence that matches what this deposit actually contains:

   > The filtered 4,002-locus by 222-individual genotype matrix, all analysis code, the complete
   > result tables and the figures are archived at https://doi.org/10.5281/zenodo.NNNNNNN and
   > developed at https://github.com/PlatypusConservation/KI-platypus. The full co-processed DArT
   > report (22,054 loci, 376 samples) is available from the corresponding author subject to the
   > permission of the owners of the mainland samples; those mainland genotypes are separately
   > available from Mijangos et al. (2022), https://doi.org/10.5281/zenodo.7039778. Raw DArTseq reads
   > are not deposited. Per-individual capture coordinates are withheld because the platypus is
   > listed Endangered in South Australia.

2. Replace the three `10.5281/zenodo.XXXXXXX` placeholders — one in `README.md` §5, one in
   `README.md` §1, one in `CITATION.cff` — and set `date-released` in `CITATION.cff` to the day you
   published. Find them with:

   ```powershell
   Select-String -Path README.md, CITATION.cff -Pattern "XXXXXXX"
   ```

3. Optional: add the Zenodo DOI badge to the top of `README.md`. Zenodo gives you the markdown on the
   record page, under **Cite as** → the badge image.

4. Commit and push the changes, and cut `v1.0.1` (§1.7) if you want the archived copy to contain its
   own DOI.

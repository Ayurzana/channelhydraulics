# channelhydraulics — Manning calculator (React + Vite)

This is a minimal React starter that implements a Manning formula calculator and is configured to deploy to GitHub Pages.

Key points
- Vite build outputs into `docs/` (configured in vite.config.js).
- base is set to `/channelhydraulics/` — update if you rename your repository.
- A GitHub Actions workflow (`.github/workflows/pages.yml`) builds and deploys automatically on push to `main`.

Local development
1. Install Node.js (>=16 recommended).
2. Clone the repo and install:
   ```
   git clone https://github.com/Ayurzana/channelhydraulics.git
   cd channelhydraulics
   npm install
   ```
3. Run the dev server:
   ```
   npm run dev
   ```
   Visit the printed localhost URL (default: http://localhost:5173).

Build & preview
- Build for production (this writes the built site into `docs/`):
  ```
  npm run build
  ```
- Preview the production build locally:
  ```
  npm run preview
  ```

Manual deployment to GitHub Pages
1. Ensure main branch has the `docs/` folder containing the production build (run `npm run build` and commit/push the `docs/` folder).
2. In GitHub: Settings → Pages → Source: choose "main branch /docs folder" and Save.
3. Wait a few minutes. Your site should be available at:
   ```
   https://Ayurzana.github.io/channelhydraulics/
   ```
   (GitHub shows the exact URL in the Pages settings.)

Automatic deployment (recommended)
- The included GitHub Actions workflow `/.github/workflows/pages.yml` builds and deploys to GitHub Pages whenever you push to `main`. No extra token is needed (it uses the repository's built-in Pages deployment permissions).

Notes and next steps
- If you change the repository name, update `base` in `vite.config.js`.
- Add more calculators as more components / routes.
- Add unit tests for calculation functions (Jest/ Vitest) and enable CI checks.
- Consider KaTeX/MathJax if you want prettier equation rendering.
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  // IMPORTANT: set base to "/<repo-name>/" if you will host on GitHub Pages for a project site
  base: '/channelhydraulics/',
  build: {
    outDir: 'docs', // build directly to `docs/` so GitHub Pages can serve from main/docs
    emptyOutDir: true
  }
});
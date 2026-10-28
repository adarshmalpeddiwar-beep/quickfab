#!/usr/bin/env python3
"""Inline src/engine.js into src/app.html and write web/index.html (a standalone static page)."""
import pathlib
root = pathlib.Path(__file__).resolve().parent.parent
app = (root / 'src/app.html').read_text()
engine = (root / 'src/engine.js').read_text()
head = '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
title_end = app.index('</title>') + len('</title>')
page = head + app[:title_end] + '\n</head>\n<body>\n' + app[title_end:].replace('<script>/*ENGINE*/</script>', '<script>\n' + engine + '\n</script>') + '\n</body>\n</html>\n'
(root / 'web/index.html').write_text(page)
print('wrote web/index.html', len(page), 'bytes')

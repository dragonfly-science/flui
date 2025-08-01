# This is what I can do
default:
  just --list

# Sync the dependencies
sync:
  uv sync

# Lint all files
lint: sync
  uv run ruff check src tests scripts

# Typecheck all files
typecheck: sync
  uv run pyright src scripts

# Run all tests
test: sync
  uv run pytest tests --benchmark-skip -ra

# Update snapshot tests
snap: sync
  uv run pytest --snapshot-update

# Save current benchmark results
benchsave: sync
  uv run pytest --benchmark-only --benchmark-autosave

# Benchmark (and compare)
bench: sync
  uv run pytest --benchmark-only --benchmark-compare

# Run all checks
check: lint typecheck test

